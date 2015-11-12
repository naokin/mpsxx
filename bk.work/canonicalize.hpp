#ifndef __MPSXX_CANONICALIZE_HPP
#define __MPSXX_CANONICALIZE_HPP

#include <legacy/QSPARSE/QSTArray.h>

namespace mpsxx {

//! Canonicalize one-site MPS
template<class Q>
void canonicalize (
        bool forward,
  const btas::QSTArray<double,3,Q>& wfnc,
        btas::QSTArray<double,3,Q>& lmps,
        btas::QSTArray<double,3,Q>& rmps,
        int M = 0)
{
  btas::STArray<double,1> s;
  btas::QSTArray<double,2,Q> g;
  if(forward) {
    btas::Gesvd<double,3,3,Q,btas::LeftArrow>(wfnc,s,lmps,g,M);
    btas::Dimm(s,g);
    btas::QSTArray<double,3,Q> temp;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,g,rmps,1.0,temp);
    rmps = temp;
  }
  else {
    btas::Gesvd<double,3,2,Q,btas::RightArrow>(wfnc,s,g,rmps,M);
    btas::Dimm(g,s);
    btas::QSTArray<double,3,Q> temp;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lmps,g,1.0,temp);
    lmps = temp;
  }
}

//! Canonicalize two-site MPS
template<class Q>
void canonicalize (
        bool forward,
  const btas::QSTArray<double,4,Q>& wfnc,
        btas::QSTArray<double,3,Q>& lmps,
        btas::QSTArray<double,3,Q>& rmps,
        int M = 0)
{
  btas::STArray<double,1> s;
  if(forward) {
    btas::Gesvd<double,4,3,Q,btas::LeftArrow>(wfnc,s,lmps,rmps,M);
    btas::Dimm(s,rmps);
  }
  else {
    btas::Gesvd<double,4,3,Q,btas::RightArrow>(wfnc,s,lmps,rmps,M);
    btas::Dimm(lmps,s);
  }
}

//! Canonicalize merged MPS
template<class Q>
void canonicalize (
        bool forward,
  const btas::QSTArray<double,2,Q>& wfnc,
        btas::QSTArray<double,2,Q>& lmat,
        btas::QSTArray<double,2,Q>& rmat,
        int M = 0)
{
  btas::STArray<double,1> s;
  if(forward) {
    btas::Gesvd<double,2,2,Q,btas::LeftArrow>(wfnc,s,lmat,rmat,M);
    btas::Dimm(s,rmat);
  }
  else {
    btas::Gesvd<double,2,2,Q,btas::RightArrow>(wfnc,s,lmat,rmat,M);
    btas::Dimm(lmat,s);
  }
}

} // namespace mpsxx

#endif // __MPSXX_CANONICALIZE_HPP
