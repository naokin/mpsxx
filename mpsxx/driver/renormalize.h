#ifndef _MPSXX_CXX11_RENORMALIZE_H
#define _MPSXX_CXX11_RENORMALIZE_H 1

#include <legacy/QSPARSE/QSDArray.h>

namespace mpsxx {

//! Renormalize operator tensor
template<class Q>
void renormalize
(bool forward, const btas::QSDArray<4, Q>& mpo0,
               const btas::QSDArray<3, Q>& opr0,
               const btas::QSDArray<3, Q>& bra0,
               const btas::QSDArray<3, Q>& ket0,
                     btas::QSDArray<3, Q>& opr1)
{
  if(forward) {
    btas::QSDArray<4, Q> scr1;
    btas::QSDcontract(1.0, opr0, btas::shape(0), bra0.conjugate(), btas::shape(0), 1.0, scr1);
    btas::QSDArray<4, Q> scr2;
    btas::QSDcontract(1.0, scr1, btas::shape(0, 2), mpo0, btas::shape(0, 1), 1.0, scr2);
    btas::QSDcontract(1.0, scr2, btas::shape(0, 2), ket0, btas::shape(0, 1), 1.0, opr1);
  }
  else {
    btas::QSDArray<4, Q> scr1;
    btas::QSDcontract(1.0, bra0.conjugate(), btas::shape(2), opr0, btas::shape(0), 1.0, scr1);
    btas::QSDArray<4, Q> scr2;
    btas::QSDcontract(1.0, scr1, btas::shape(1, 2), mpo0, btas::shape(1, 3), 1.0, scr2);
    btas::QSDcontract(1.0, scr2, btas::shape(3, 1), ket0, btas::shape(1, 2), 1.0, opr1);
  }
}

//! Renormalize overlap matrix
template<class Q>
void renormalize
(bool forward, const btas::QSDArray<2, Q>& sov0,
               const btas::QSDArray<3, Q>& bra0,
               const btas::QSDArray<3, Q>& ket0,
                     btas::QSDArray<2, Q>& sov1)
{
  if(forward) {
    btas::QSDArray<3, Q> scr1;
    btas::QSDgemm(btas::Trans, btas::NoTrans, 1.0, sov0, bra0.conjugate(), 1.0, scr1);
    btas::QSDgemm(btas::Trans, btas::NoTrans, 1.0, scr1, ket0,             1.0, sov1);
  }
  else {
    btas::QSDArray<3, Q> scr1;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, bra0.conjugate(), sov0, 1.0, scr1);
    btas::QSDgemm(btas::NoTrans, btas::Trans,   1.0, scr1, ket0,             1.0, sov1);
  }
}

//! Renormalize from merged LS-block
template<class Q>
void renormalize
(bool forward, const btas::QSDArray<3, Q>& opr0,
               const btas::QSDArray<2, Q>& bra0,
               const btas::QSDArray<2, Q>& ket0,
                     btas::QSDArray<3, Q>& opr1)
{
  if(forward) {
    btas::QSDArray<3, Q> scr1;
    btas::QSDgemm(btas::Trans,   btas::NoTrans, 1.0, bra0.conjugate(), opr0, 1.0, scr1);
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, scr1, ket0,             1.0, opr1);
  }
  else {
    btas::QSDArray<3, Q> scr1;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, bra0.conjugate(), opr0, 1.0, scr1);
    btas::QSDgemm(btas::NoTrans, btas::Trans,   1.0, scr1, ket0,             1.0, opr1);
  }
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_RENORMALIZE_H
