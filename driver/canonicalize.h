#ifndef _MPSXX_CXX11_CANONICALIZE_H
#define _MPSXX_CXX11_CANONICALIZE_H 1

#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx {

//! Canonicalize one-site MPS
template<class Q>
void canonicalize
(bool forward, const btas::QSDArray<3, Q>& wfn0,
                     btas::QSDArray<3, Q>& mps0, int M = 0)
{
  if(forward) {
    btas:: SDArray<1>    s;
    btas::QSDArray<2, Q> v;
    btas::QSDgesvd(btas::LeftArrow,  wfn0, s, mps0, v, M);
  }
  else {
    btas:: SDArray<1>    s;
    btas::QSDArray<2, Q> u;
    btas::QSDgesvd(btas::RightArrow, wfn0, s, u, mps0, M);
  }
}

//! Canonicalize two-site MPS
template<class Q>
void canonicalize
(bool forward, const btas::QSDArray<4, Q>& wfnx,
                     btas::QSDArray<3, Q>& mps0,
                     btas::QSDArray<3, Q>& wfn1, int M = 0)
{
  if(forward) {
    btas::SDArray<1> s;
    btas::QSDgesvd(btas::LeftArrow,  wfnx, s, mps0, wfn1, M);
    btas::SDdidm(s, wfn1);
  }
  else {
    btas::SDArray<1> s;
    btas::QSDgesvd(btas::RightArrow, wfnx, s, wfn1, mps0, M);
    btas::SDdimd(wfn1, s);
  }
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_CANONICALIZE_H
