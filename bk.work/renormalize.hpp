#ifndef __MPSXX_RENORMALIZE_HPP
#define __MPSXX_RENORMALIZE_HPP

#include <legacy/QSPARSE/QSTArray.h>

namespace mpsxx {

//! Renormalize operator tensor
template<class Q>
void renormalize (
        bool forward,
  const btas::QSTArray<double,4,Q>& mpo0,
  const btas::QSTArray<double,3,Q>& opr0,
  const btas::QSTArray<double,3,Q>& bra0,
  const btas::QSTArray<double,3,Q>& ket0,
        btas::QSTArray<double,3,Q>& opr1)
{
  if(forward) {
    btas::QSTArray<double,4,Q> scr1;
    btas::Contract(1.0,opr0,btas::shape(0),bra0.conjugate(),btas::shape(0),1.0,scr1);
    btas::QSTArray<double,4,Q> scr2;
    btas::Contract(1.0,scr1,btas::shape(0,2),mpo0,btas::shape(0,1),1.0,scr2);
    btas::Contract(1.0,scr2,btas::shape(0,2),ket0,btas::shape(0,1),1.0,opr1);
  }
  else {
    btas::QSTArray<double,4,Q> scr1;
    btas::Contract(1.0,bra0.conjugate(),btas::shape(2),opr0,btas::shape(0),1.0,scr1);
    btas::QSTArray<double,4,Q> scr2;
    btas::Contract(1.0,scr1,btas::shape(1,2),mpo0,btas::shape(1,3),1.0,scr2);
    btas::Contract(1.0,scr2,btas::shape(3,1),ket0,btas::shape(1,2),1.0,opr1);
  }
}

//! Renormalize overlap matrix
template<class Q>
void renormalize (
        bool forward,
  const btas::QSTArray<double,2,Q>& sov0,
  const btas::QSTArray<double,3,Q>& bra0,
  const btas::QSTArray<double,3,Q>& ket0,
        btas::QSTArray<double,2,Q>& sov1)
{
  if(forward) {
    btas::QSTArray<double,3,Q> scr1;
    btas::Gemm(btas::Trans,btas::NoTrans,1.0,sov0,bra0.conjugate(),1.0,scr1);
    btas::Gemm(btas::Trans,btas::NoTrans,1.0,scr1,ket0,1.0,sov1);
  }
  else {
    btas::QSTArray<double,3,Q> scr1;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,bra0.conjugate(),sov0,1.0,scr1);
    btas::Gemm(btas::NoTrans,btas::Trans,1.0,scr1,ket0,1.0,sov1);
  }
}

//! Renormalize from merged LS-block
template<class Q>
void renormalize (
        bool forward,
  const btas::QSTArray<double,3,Q>& opr0,
  const btas::QSTArray<double,2,Q>& bra0,
  const btas::QSTArray<double,2,Q>& ket0,
        btas::QSTArray<double,3,Q>& opr1)
{
  if(forward) {
    btas::QSTArray<double,3,Q> scr1;
    btas::Gemm(btas::Trans,btas::NoTrans,1.0,bra0.conjugate(),opr0,1.0,scr1);
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,scr1,ket0,1.0,opr1);
  }
  else {
    btas::QSTArray<double,3,Q> scr1;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,bra0.conjugate(),opr0,1.0,scr1);
    btas::Gemm(btas::NoTrans,btas::Trans,1.0,scr1,ket0,1.0,opr1);
  }
}

} // namespace mpsxx

#endif // __MPSXX_RENORMALIZE_HPP
