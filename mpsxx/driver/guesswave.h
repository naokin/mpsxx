#ifndef _MPSXX_CXX11_COMPUTE_GUESSWAVE_H
#define _MPSXX_CXX11_COMPUTE_GUESSWAVE_H 1

#include <legacy/QSPARSE/QSDArray.h>

namespace mpsxx {

//! Compute guess wavefunction for next site
/*! The gauge of matrix product is moved */
template<class Q>
void compute_guesswave
(bool forward, const btas::QSDArray<3, Q>& mps0,
               const btas::QSDArray<3, Q>& wfn0,
               const btas::QSDArray<3, Q>& mps1,
                     btas::QSDArray<3, Q>& wfn1)
{
  if(forward) {
    btas::QSDArray<2, Q> lres;
    btas::QSDgemm(btas::ConjTrans, btas::NoTrans, 1.0, mps0, wfn0, 1.0, lres);
    btas::QSDgemm(  btas::NoTrans, btas::NoTrans, 1.0, lres, mps1, 1.0, wfn1);
  }
  else {
    btas::QSDArray<2, Q> rres;
    btas::QSDgemm(btas::NoTrans, btas::ConjTrans, 1.0, wfn0, mps0, 1.0, rres);
    btas::QSDgemm(btas::NoTrans,   btas::NoTrans, 1.0, mps1, rres, 1.0, wfn1);
  }
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_COMPUTE_GUESSWAVE_H
