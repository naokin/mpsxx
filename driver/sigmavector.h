#ifndef _MPSXX_CXX11_SIGMAVECTOR_H
#define _MPSXX_CXX11_SIGMAVECTOR_H 1

#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx {

//! Compute sigmavector: | phi > = H | psi >
template<class Q>
void compute_sigmavector
(const btas::QSDArray<4, Q>& mpo0,
 const btas::QSDArray<3, Q>& lopr,
 const btas::QSDArray<3, Q>& ropr,
 const btas::QSDArray<3, Q>& wfn0,
       btas::QSDArray<3, Q>& sgv0)
{
  btas::QSDArray<4, Q> scr1;
  btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, lopr, wfn0, 1.0, scr1);
  btas::QSDArray<4, Q> scr2;
  btas::QSDcontract(1.0, scr1, btas::shape(1, 2), mpo0, btas::shape(0, 2), 1.0, scr2);
  btas::QSDcontract(1.0, scr2, btas::shape(3, 1), ropr, btas::shape(1, 2), 1.0, sgv0);
}

//! Compute sigmavector: < phi | = < psi | H^*
template<class Q>
void compute_sigmavector_conj
(const btas::QSDArray<4>& mpo0,
 const btas::QSDArray<3>& lopr,
 const btas::QSDArray<3>& ropr,
 const btas::QSDArray<3>& wfn0,
       btas::QSDArray<3>& sgv0)
{
  btas::QSDArray<4, Q> scr1;
  btas::QSDgemm(btas::Trans, btas::ConjTrans, 1.0, lopr, wfn0, 1.0, scr1);
  btas::QSDArray<4, Q> scr2;
  btas::QSDcontract(1.0, scr1, btas::shape(0, 2), mpo0, btas::shape(0, 1), 1.0, scr2);
  btas::QSDcontract(1.0, scr2, btas::shape(1, 3), ropr, btas::shape(0, 1), 1.0, sgv0.conjugate());
}

//! Compute sigmavector for two-site algorithm: | phi > = H | psi >
template<class Q>
void compute_sigmavector
(const btas::QSDArray<4, Q>& lmpo,
 const btas::QSDArray<4, Q>& rmpo,
 const btas::QSDArray<3, Q>& lopr,
 const btas::QSDArray<3, Q>& ropr,
 const btas::QSDArray<4, Q>& wfn0,
       btas::QSDArray<4, Q>& sgv0)
{
  btas::QSDArray<5, Q> scr1;
  btas::QSDcontract(1.0, lopr, btas::shape(2),    wfn0, btas::shape(0),    1.0, scr1);
  btas::QSDArray<5, Q> scr2;
  btas::QSDcontract(1.0, scr1, btas::shape(1, 2), lmpo, btas::shape(0, 2), 1.0, scr2);
  btas::QSDArray<5, Q> scr3;
  btas::QSDcontract(1.0, scr2, btas::shape(4, 1), rmpo, btas::shape(0, 2), 1.0, scr3);
  btas::QSDcontract(1.0, scr3, btas::shape(4, 1), ropr, btas::shape(1, 2), 1.0, sgv0);
}

//! Compute sigmavector for merged algorithm: | phi > = H | psi >
template<class Q>
void compute_sigmavector
(const btas::QSDArray<3, Q>& lopr,
 const btas::QSDArray<3, Q>& ropr,
 const btas::QSDArray<2, Q>& wfn0,
       btas::QSDArray<2, Q>& sgv0)
{
  btas::QSDArray<3, Q> scr1;
  btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, lopr, wfn0, 1.0, scr1);
  btas::QSDgemm(btas::NoTrans, btas::Trans,   1.0, scr1, ropr, 1.0, sgv0);
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_SIGMAVECTOR_H
