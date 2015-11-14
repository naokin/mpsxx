#ifndef __MPSXX_COMPUTE_GUESSWAVE_HPP
#define __MPSXX_COMPUTE_GUESSWAVE_HPP

#include <btas/QSPARSE/QSTArray.h>

namespace mpsxx {

/// Compute guess wavefunction for next site
/// The gauge of matrix product is tranferred.
template<class Q>
void compute_guesswave (
        bool forward,
  const btas::QSTArray<double,3,Q>& mps0,
  const btas::QSTArray<double,3,Q>& wfn0,
  const btas::QSTArray<double,3,Q>& mps1,
        btas::QSTArray<double,3,Q>& wfn1)
{
  if(forward) {
    // mps0 : L(i)
    // wfn0 : C(i)
    // mps1 : R(i+1)
    // wfn1 : C(i+1)
    // C(i+1) = L(i)^T * C(i) * R(i+1)
    btas::QSTArray<double,2,Q> gmat;
    btas::Gemm(btas::ConjTrans, btas::NoTrans,1.0,mps0,wfn0,1.0,gmat);
    btas::Gemm(  btas::NoTrans, btas::NoTrans,1.0,gmat,mps1,1.0,wfn1);
  }
  else {
    // mps0 : R(i+1)
    // wfn0 : C(i+1)
    // mps1 : L(i)
    // wfn1 : C(i)
    // C(i) = L(i) * C(i+1) * R(i+1)^T
    btas::QSTArray<double,2,Q> gmat;
    btas::Gemm(btas::NoTrans, btas::ConjTrans,1.0,wfn0,mps0,1.0,gmat);
    btas::Gemm(btas::NoTrans,   btas::NoTrans,1.0,mps1,gmat,1.0,wfn1);
  }
}

} // namespace mpsxx

#endif // __MPSXX_COMPUTE_GUESSWAVE_HPP
