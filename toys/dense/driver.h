#ifndef _PROTOTYPE_DRIVER_H
#define _PROTOTYPE_DRIVER_H 1

#include <btas/DArray.h>
#include <btas/Dlapack.h>

namespace prototype
{

template<int NA, int NU>
void Decompose
(const btas::DArray<NA>& a, btas::DArray<1>& s, btas::DArray<NU>& u, btas::DArray<NA-NU+2>& vt, int M = 0)
{
  btas::DArray<1>        s_tmp;
  btas::DArray<NU>       u_tmp;
  btas::DArray<NA-NU+2> vt_tmp;
  btas::Dgesvd(a, s_tmp, u_tmp, vt_tmp);
  if(M == 0 || M >= s_tmp.size()) {
    btas::Dcopy( s_tmp,  s);
    btas::Dcopy( u_tmp,  u);
    btas::Dcopy(vt_tmp, vt);
    return;
  }
  btas::TinyVector<int, NU>       u_shape( u_tmp.shape());
  btas::TinyVector<int, NA-NU+2> vt_shape(vt_tmp.shape());
   u_shape[NU-1] = M;
  vt_shape[0]    = M;

   s.resize(M);
   u.resize(u_shape);
  vt.resize(vt_shape);

  for(int i = 0; i < M; ++i)
    s(i) = s_tmp(i);

  for(int i = 0; i < NU;      ++i) -- u_shape[i];
  for(int i = 0; i < NA-NU+2; ++i) --vt_shape[i];

  btas::TinyVector<int, NU>       u_lbound(0);
  btas::TinyVector<int, NA-NU+2> vt_lbound(0);

  btas::Dcopy_direct( u_tmp(btas::RectDomain<NU>     ( u_lbound,  u_shape)),  u);
  btas::Dcopy_direct(vt_tmp(btas::RectDomain<NA-NU+2>(vt_lbound, vt_shape)), vt);
}

void ComputeGuess
(bool forward, const btas::DArray<3>& mps0,
               const btas::DArray<3>& wfn0,
               const btas::DArray<3>& mps1,
                     btas::DArray<3>& wfn1);

void Canonicalize
(bool forward, const btas::DArray<3>& wfn0,
                     btas::DArray<3>& mps0, int M = 0);

void Canonicalize
(bool forward, const btas::DArray<3>& wfn0,
                     btas::DArray<3>& mps0,
                     btas::DArray<3>& mps1);

void Canonicalize
(bool forward, const btas::DArray<4>& wfnx,
                     btas::DArray<3>& mps0,
                     btas::DArray<3>& wfn1, int M = 0);

void Renormalize
(bool forward, const btas::DArray<4>& mpo0,
               const btas::DArray<3>& opr0,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket0,
                     btas::DArray<3>& opr1);

void Renormalize
(bool forward, const btas::DArray<2>& opr0,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket0,
                     btas::DArray<2>& opr1);

void ComputeInverseGauge
(              const btas::DArray<2>& gauge_0,
                     btas::DArray<2>& gauge_i);

void Normalize
(                    btas::DArray<3>& wfn0);

void Orthogonalize
(              const btas::DArray<3>& wfn0,
                     btas::DArray<3>& wfn1);

void Orthogonalize
(bool forward, const btas::DArray<3>& mps0,
                     btas::DArray<3>& wfn1);


void ComputeDiagonal
(              const btas::DArray<4>& mpo0,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
                     btas::DArray<3>& diag);

void ComputeDiagonal
(              const btas::DArray<4>& lmpo,
               const btas::DArray<4>& rmpo,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
                     btas::DArray<4>& diag);

void ComputeSigmaVector
(              const btas::DArray<4>& mpo0,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
               const btas::DArray<3>& wfn0,
                     btas::DArray<3>& sgv0);

void ComputeSigmaVectorConj
(              const btas::DArray<4>& mpo0,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
               const btas::DArray<3>& wfn0,
                     btas::DArray<3>& sgv0);

void ComputeSigmaVector
(              const btas::DArray<4>& lmpo,
               const btas::DArray<4>& rmpo,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
               const btas::DArray<4>& wfn0,
                     btas::DArray<4>& sgv0);

int ComputeOrthogonalTransform
(              const btas::DArray<2>& s_matrix,
                     btas::DArray<2>& u_matrix);

void TransformedMatrix
(              const btas::DArray<2>& h_matrix,
               const btas::DArray<2>& u_matrix,
                     btas::DArray<2>& h_transm);

int ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& s_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alpha);

int ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm);

namespace TEST {

int ComputeOrthogonalTransform
(              const btas::DArray<2>& s_matrix,
                     btas::DArray<2>& u_matrix);

int ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm);

};

};

#endif // _PROTOTYPE_DRIVER_H
