#ifndef __MPSXX_DMRG_DENSE_DRIVER_HPP
#define __MPSXX_DMRG_DENSE_DRIVER_HPP

#include "MPX.h"

namespace mpsxx {

template<typename T>
void compute_guess_wave (bool forward, const MPS<T>& mps0, const MPS<T>& wfn0, const MPS<T>& mps1, MPS<T>& wfn1)
{
  if(forward) {
    btas::TArray<T,2> gmat;
    btas::Gemm(CblasTrans,  CblasNoTrans,1.0,mps0,wfn0,1.0,gmat);
    btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,gmat,mps1,1.0,wfn1);
  }
  else {
    btas::TArray<T,2> gmat;
    btas::Gemm(CblasNoTrans,  CblasTrans,1.0,wfn0,mps0,1.0,gmat);
    btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,mps1,gmat,1.0,wfn1);
  }
}

template<typename T, size_t N>
void canonicalize (bool forward, const btas::TArray<T,N>& wfn0, btas::TArray<T,3>& mps0, btas::TArray<T,N-1>& res0, size_t M = 0)
{
  const size_t K = N-1;

  btas::TArray<typename btas::remove_complex<T>::type,1> s;
  btas::TArray<T,3> mpst;
  btas::TArray<T,K> rest;

  if(forward)
    btas::Gesvd('S','S',wfn0,s,mpst,rest);
  else
    btas::Gesvd('S','S',wfn0,s,rest,mpst);

  size_t M0 = s.size();

  if(M == 0) {
    for(; M < M0; ++M) if(s(M) < 1.0e-12) break;
  }
  if(M > M0) {
    M = M0;
  }

  double factor = 0.0;
  for(size_t i = 0; i < M; ++i) factor += s(i)*s(i);

  factor = 1.0/sqrt(factor);

//if(M < M0) {
    if(forward) {
      auto Dmps = mpst.shape();
      Dmps[2] = M;
      size_t L = Dmps[0]*Dmps[1];

      mps0.resize(Dmps);
      {
        const T* x = mpst.data();
              T* y = mps0.data();
        for(size_t i = 0; i < L; ++i)
          for(size_t j = 0; j < M; ++j) y[i*M+j] = x[i*M0+j];
      }

      auto Dres = rest.shape();
      Dres[0] = M;
      size_t R = 1; for(size_t i = 1; i < K; ++i) R *= Dres[i];

      res0.resize(Dres);
      {
        const T* x = rest.data();
              T* y = res0.data();
        for(size_t i = 0; i < M; ++i)
          for(size_t j = 0; j < R; ++j) y[i*R+j] = factor*s(i)*x[i*R+j];
      }
    }
    else {
      auto Dres = rest.shape();
      Dres[K-1] = M;
      size_t L = 1; for(size_t i = 0; i < K-1; ++i) L *= Dres[i];

      res0.resize(Dres);
      {
        const T* x = rest.data();
              T* y = res0.data();
        for(size_t i = 0; i < L; ++i)
          for(size_t j = 0; j < M; ++j) y[i*M+j] = x[i*M0+j]*factor*s(j);
      }

      auto Dmps = mpst.shape();
      Dmps[0] = M;
      size_t R = Dmps[1]*Dmps[2];

      mps0.resize(Dmps);
      {
        const T* x = mpst.data();
              T* y = mps0.data();
        for(size_t i = 0; i < M; ++i)
          for(size_t j = 0; j < R; ++j) y[i*R+j] = x[i*R+j];
      }
    }
//}
//else {
//  mps0 = mpst;
//  res0 = rest;
//}
}

template<typename T>
void renormalize (bool forward, const MPO<T>& mpo0, const BLOCK<T>& opr0, const MPS<T>& bra0, const MPS<T>& ket0, BLOCK<T>& opr1)
{
  using btas::shape;

  if(forward) {
    btas::TArray<T,4> tmp1;
    btas::Contract(1.0,opr0,shape(0),  bra0,shape(0),  1.0,tmp1);
    btas::TArray<T,4> tmp2;
    btas::Contract(1.0,tmp1,shape(0,2),mpo0,shape(0,1),1.0,tmp2);
    btas::Contract(1.0,tmp2,shape(0,2),ket0,shape(0,1),1.0,opr1);
  }
  else {
    btas::TArray<T,4> tmp1;
    btas::Contract(1.0,bra0,shape(2),  opr0,shape(0),  1.0,tmp1);
    btas::TArray<T,4> tmp2;
    btas::Contract(1.0,tmp1,shape(1,2),mpo0,shape(1,3),1.0,tmp2);
    btas::Contract(1.0,tmp2,shape(3,1),ket0,shape(1,2),1.0,opr1);
  }
}

template<typename T>
void compute_diagonal_elements (const MPO<T>& mpo0, const BLOCK<T>& lopr, const BLOCK<T>& ropr, MPS<T>& diag)
{
  using btas::shape;

  btas::TArray<T,3> mpo0_diag(mpo0.shape(0),mpo0.shape(1),mpo0.shape(3));
  btas::TArray<T,2> lopr_diag(lopr.shape(0),lopr.shape(1));
  btas::TArray<T,2> ropr_diag(ropr.shape(0),ropr.shape(1));

  for(size_t i = 0; i < mpo0.shape(0); ++i)
    for(size_t j = 0; j < mpo0.shape(1); ++j)
      for(size_t k = 0; k < mpo0.shape(3); ++k)
        mpo0_diag(i,j,k) = mpo0(i,j,j,k);

  for(size_t i = 0; i < lopr.shape(0); ++i)
    for(size_t j = 0; j < lopr.shape(1); ++j)
      lopr_diag(i,j) = lopr(i,j,i);

  for(size_t i = 0; i < ropr.shape(0); ++i)
    for(size_t j = 0; j < ropr.shape(1); ++j)
      ropr_diag(i,j) = ropr(i,j,i);

  btas::TArray<T,3> tmp1;
  btas::Contract(1.0,lopr_diag,shape(1),mpo0_diag,shape(0),1.0,tmp1);
  btas::Contract(1.0,tmp1,     shape(2),ropr_diag,shape(1),1.0,diag);
}

template<typename T>
void compute_diagonal_elements (const MPO<T>& lmpo, const MPO<T>& rmpo, const BLOCK<T>& lopr, const BLOCK<T>& ropr, MPS<T>& diag)
{
  using btas::shape;

  btas::TArray<T,3> lmpo_diag(lmpo.shape(0),lmpo.shape(1),lmpo.shape(3));
  btas::TArray<T,3> rmpo_diag(rmpo.shape(0),rmpo.shape(1),rmpo.shape(3));
  btas::TArray<T,2> lopr_diag(lopr.shape(0),lopr.shape(1));
  btas::TArray<T,2> ropr_diag(ropr.shape(0),ropr.shape(1));

  for(size_t i = 0; i < lmpo.shape(0); ++i)
    for(size_t j = 0; j < lmpo.shape(1); ++j)
      for(size_t k = 0; k < lmpo.shape(3); ++k)
        lmpo_diag(i,j,k) = lmpo(i,j,j,k);

  for(size_t i = 0; i < rmpo.shape(0); ++i)
    for(size_t j = 0; j < rmpo.shape(1); ++j)
      for(size_t k = 0; k < rmpo.shape(3); ++k)
        rmpo_diag(i,j,k) = rmpo(i,j,j,k);

  for(size_t i = 0; i < lopr.shape(0); ++i)
    for(size_t j = 0; j < lopr.shape(1); ++j)
      lopr_diag(i,j) = lopr(i,j,i);

  for(size_t i = 0; i < ropr.shape(0); ++i)
    for(size_t j = 0; j < ropr.shape(1); ++j)
      ropr_diag(i,j) = ropr(i,j,i);

  btas::TArray<T,3> tmp1;
  btas::Contract(1.0,lopr_diag,shape(1),lmpo_diag,shape(0),1.0,tmp1);
  btas::TArray<T,4> tmp2;
  btas::Contract(1.0,tmp1,     shape(2),rmpo_diag,shape(0),1.0,tmp2);
  btas::Contract(1.0,tmp2,     shape(3),ropr_diag,shape(1),1.0,diag);
}

/// compute H | psi >
template<typename T>
void compute_sigma_vector (const MPO<T>& mpo0, const BLOCK<T>& lopr, const BLOCK<T>& ropr, const MPS<T>& wfn0, MPS<T>& sgv0)
{
  using btas::shape;

  btas::TArray<T,4> tmp1;
  btas::Contract(1.0,lopr,shape(2),  wfn0,shape(0),  1.0,tmp1);
  btas::TArray<T,4> tmp2;
  btas::Contract(1.0,tmp1,shape(1,2),mpo0,shape(0,2),1.0,tmp2);
  btas::Contract(1.0,tmp2,shape(3,1),ropr,shape(1,2),1.0,sgv0);
}

template<typename T>
void compute_sigma_vector (const MPO<T>& lmpo, const MPO<T>& rmpo, const BLOCK<T>& lopr, const BLOCK<T>& ropr, const MPS<T>& wfn0, MPS<T>& sgv0)
{
  using btas::shape;

  btas::TArray<T,5> tmp1;
  btas::Contract(1.0,lopr,shape(2),  wfn0,shape(0),  1.0,tmp1);
  btas::TArray<T,5> tmp2;
  btas::Contract(1.0,tmp1,shape(1,2),lmpo,shape(0,2),1.0,tmp2);
  btas::TArray<T,5> tmp3;
  btas::Contract(1.0,tmp2,shape(4,1),rmpo,shape(0,2),1.0,tmp3);
  btas::Contract(1.0,tmp3,shape(4,1),ropr,shape(1,2),1.0,sgv0);
}

} // namespace mpsxx

#endif // __MPSXX_DMRG_DENSE_DRIVER_HPP
