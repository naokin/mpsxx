#ifndef __MPSXX_DRIVER_CANONICALIZE_HPP
#define __MPSXX_DRIVER_CANONICALIZE_HPP

#include <legacy/QSPARSE/QSDArray.h>

namespace mpsxx {

/// Canonicalize one-site MPS
template<class Q> void canonicalize (
        bool forward,
  const btas::QSDArray<3, Q>& wfn0,
        btas::QSDArray<3, Q>& mps0,
        btas::QSDArray<2, Q>& wfn1, int M = 0)
{
  btas::SDArray<1> s;
  if(forward) {
    btas::QSDgesvd(btas::LeftArrow,  wfn0, s, mps0, wfn1, M);
    btas::Dimm(s, wfn1);
  }
  else {
    btas::QSDgesvd(btas::RightArrow, wfn0, s, wfn1, mps0, M);
    btas::Dimm(wfn1, s);
  }
}

/// Canonicalize two-site MPS
template<class Q> void canonicalize (
        bool forward,
  const btas::QSDArray<4, Q>& wfnx,
        btas::QSDArray<3, Q>& mps0,
        btas::QSDArray<3, Q>& wfn1, int M = 0)
{
  btas::SDArray<1> s;
  if(forward) {
    btas::QSDgesvd(btas::LeftArrow,  wfnx, s, mps0, wfn1, M);
    btas::Dimm(s, wfn1);
  }
  else {
    btas::QSDgesvd(btas::RightArrow, wfnx, s, wfn1, mps0, M);
    btas::Dimm(wfn1, s);
  }
}

/// Canonicalize merged MPS
template<class Q> void canonicalize (
        bool forward,
  const btas::QSDArray<2, Q>& wfnx,
        btas::QSDArray<2, Q>& mps0,
        btas::QSDArray<2, Q>& wfn1, int M = 0)
{
  btas::SDArray<1> s;
  if(forward) {
    btas::QSDgesvd(btas::LeftArrow,  wfnx, s, mps0, wfn1, M);
    btas::Dimm(s, wfn1);
  }
  else {
    btas::QSDgesvd(btas::RightArrow, wfnx, s, wfn1, mps0, M);
    btas::Dimm(wfn1, s);
  }
}

} // namespace mpsxx

#endif // __MPSXX_DRIVER_CANONICALIZE_HPP
