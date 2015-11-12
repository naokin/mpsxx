#ifndef __MPSXX_COMPUTE_SIGMA_VECTOR_HPP
#define __MPSXX_COMPUTE_SIGMA_VECTOR_HPP

#include <legacy/QSPARSE/QSTArray.h>

namespace mpsxx {

//! Compute sigmavector: | phi > = H | psi >
template<class Q>
void compute_sigma_vector (
  const btas::QSTArray<double,4,Q>& mpo0,
  const btas::QSTArray<double,3,Q>& lopr,
  const btas::QSTArray<double,3,Q>& ropr,
  const btas::QSTArray<double,3,Q>& wfn0,
        btas::QSTArray<double,3,Q>& sgv0)
{
  btas::QSTArray<double,4,Q> scr1;
  btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lopr,wfn0,1.0,scr1);
  btas::QSTArray<double,4,Q> scr2;
  btas::Contract(1.0,scr1,btas::shape(1,2),mpo0,btas::shape(0,2),1.0,scr2);
  btas::Contract(1.0,scr2,btas::shape(3,1),ropr,btas::shape(1,2),1.0,sgv0);
}

//! Compute sigmavector: | phi > = S | psi >
template<class Q>
void compute_sigma_vector (
  const btas::QSTArray<double,2,Q>& lopr,
  const btas::QSTArray<double,2,Q>& ropr,
  const btas::QSTArray<double,3,Q>& wfn0,
        btas::QSTArray<double,3,Q>& sgv0)
{
  btas::QSTArray<double,3,Q> scr1;
  btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lopr,wfn0,1.0,scr1);
  btas::Gemm(btas::NoTrans,btas::Trans,1.0,scr1,ropr,1.0,sgv0);
}

//! Compute sigmavector: < phi | = < psi | H^*
template<class Q>
void compute_sigma_vector_conj (
  const btas::QSTArray<double,4,Q>& mpo0,
  const btas::QSTArray<double,3,Q>& lopr,
  const btas::QSTArray<double,3,Q>& ropr,
  const btas::QSTArray<double,3,Q>& wfn0,
        btas::QSTArray<double,3,Q>& sgv0)
{
  btas::QSTArray<double,4,Q> scr1;
  btas::Gemm(btas::Trans,btas::ConjTrans,1.0,lopr,wfn0,1.0,scr1);
  btas::QSTArray<double,4,Q> scr2;
  btas::Contract(1.0,scr1,btas::shape(0,2),mpo0,btas::shape(0,1),1.0,scr2);
  btas::Contract(1.0,scr2,btas::shape(1,3),ropr,btas::shape(0,1),1.0,sgv0.conjugate());
}

//! Compute sigmavector for two-site algorithm: | phi > = H | psi >
template<class Q>
void compute_sigma_vector (
  const btas::QSTArray<double,4,Q>& lmpo,
  const btas::QSTArray<double,4,Q>& rmpo,
  const btas::QSTArray<double,3,Q>& lopr,
  const btas::QSTArray<double,3,Q>& ropr,
  const btas::QSTArray<double,4,Q>& wfn0,
        btas::QSTArray<double,4,Q>& sgv0)
{
  btas::QSTArray<double,5,Q> scr1;
  btas::Contract(1.0,lopr,btas::shape(2),wfn0,btas::shape(0),1.0,scr1);
  btas::QSTArray<double,5,Q> scr2;
  btas::Contract(1.0,scr1,btas::shape(1,2),lmpo,btas::shape(0,2),1.0,scr2);
  btas::QSTArray<double,5,Q> scr3;
  btas::Contract(1.0,scr2,btas::shape(4,1),rmpo,btas::shape(0,2),1.0,scr3);
  btas::Contract(1.0,scr3,btas::shape(4,1),ropr,btas::shape(1,2),1.0,sgv0);
}

//! Compute sigmavector: | phi > = S | psi >
template<class Q>
void compute_sigma_vector (
  const btas::QSTArray<double,2,Q>& lopr,
  const btas::QSTArray<double,2,Q>& ropr,
  const btas::QSTArray<double,4,Q>& wfn0,
        btas::QSTArray<double,4,Q>& sgv0)
{
  btas::QSTArray<double,4,Q> scr1;
  btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lopr,wfn0,1.0,scr1);
  btas::Gemm(btas::NoTrans,btas::Trans,1.0,scr1,ropr,1.0,sgv0);
}

//! Compute sigmavector for merged algorithm: | phi > = H | psi >
template<class Q>
void compute_sigma_vector (
  const btas::QSTArray<double,3,Q>& lopr,
  const btas::QSTArray<double,3,Q>& ropr,
  const btas::QSTArray<double,2,Q>& wfn0,
        btas::QSTArray<double,2,Q>& sgv0)
{
  btas::QSTArray<double,3,Q> scr1;
  btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lopr,wfn0,1.0,scr1);
  btas::Gemm(btas::NoTrans,btas::Trans,1.0,scr1,ropr,1.0,sgv0);
}

//! Compute sigmavector: | phi > = S | psi >
template<class Q>
void compute_sigma_vector (
  const btas::QSTArray<double,2,Q>& lopr,
  const btas::QSTArray<double,2,Q>& ropr,
  const btas::QSTArray<double,2,Q>& wfn0,
        btas::QSTArray<double,2,Q>& sgv0)
{
  btas::QSTArray<double,2,Q> scr1;
  btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lopr,wfn0,1.0,scr1);
  btas::Gemm(btas::NoTrans,btas::Trans,1.0,scr1,ropr,1.0,sgv0);
}

} // namespace mpsxx

#endif // __MPSXX_COMPUTE_SIGMA_VECTOR_HPP
