#ifndef __MPSXX_COMPUTE_DIAGONAL_ELEMENTS_H
#define __MPSXX_COMPUTE_DIAGONAL_ELEMENTS_H

#include <legacy/QSPARSE/QSTArray.h>

namespace mpsxx {

//! Compute diagonal H elements for one-site algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSTArray<double,4, Q>& mpo0,
 const btas::QSTArray<double,3, Q>& lopr,
 const btas::QSTArray<double,3, Q>& ropr,
       btas::QSTArray<double,3, Q>& diag)
{
  btas::STArray<double,3> mpo0_diag;
  btas::STArray<double,2> lopr_diag;
  btas::STArray<double,2> ropr_diag;

  btas::Tie(mpo0, btas::shape(1, 2), mpo0_diag);
  btas::Tie(lopr, btas::shape(0, 2), lopr_diag);
  btas::Tie(ropr, btas::shape(0, 2), ropr_diag);

  btas::STArray<double,3> scr1;
  btas::Contract(1.0, lopr_diag, btas::shape(1), mpo0_diag, btas::shape(0), 1.0, scr1);
  btas::STArray<double,3> scr2;
  btas::Contract(1.0, scr1,      btas::shape(2), ropr_diag, btas::shape(1), 1.0, scr2);
  btas::Copy(scr2, diag, 1); // preserve quantum number of diag
}

//! Compute diagonal H elements for two-site algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSTArray<double,4, Q>& lmpo,
 const btas::QSTArray<double,4, Q>& rmpo,
 const btas::QSTArray<double,3, Q>& lopr,
 const btas::QSTArray<double,3, Q>& ropr,
       btas::QSTArray<double,4, Q>& diag)
{
  btas::STArray<double,3> lmpo_diag;
  btas::STArray<double,3> rmpo_diag;
  btas::STArray<double,2> lopr_diag;
  btas::STArray<double,2> ropr_diag;

  btas::Tie(lmpo, btas::shape(1, 2), lmpo_diag);
  btas::Tie(rmpo, btas::shape(1, 2), rmpo_diag);
  btas::Tie(lopr, btas::shape(0, 2), lopr_diag);
  btas::Tie(ropr, btas::shape(0, 2), ropr_diag);

  btas::STArray<double,3> scr1;
  btas::Contract(1.0, lopr_diag, btas::shape(1), lmpo_diag, btas::shape(0), 1.0, scr1);
  btas::STArray<double,4> scr2;
  btas::Contract(1.0, scr1,      btas::shape(2), rmpo_diag, btas::shape(0), 1.0, scr2);
  btas::STArray<double,4> scr3;
  btas::Contract(1.0, scr2,      btas::shape(3), ropr_diag, btas::shape(1), 1.0, scr3);
  btas::Copy(scr3, diag, 1); // preserve quantum number of diag
}

//! Compute diagonal H elements for merged algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSTArray<double,3, Q>& lopr,
 const btas::QSTArray<double,3, Q>& ropr,
       btas::QSTArray<double,2, Q>& diag)
{
  btas::STArray<double,2> lopr_diag;
  btas::STArray<double,2> ropr_diag;

  btas::Tie(lopr, btas::shape(0, 2), lopr_diag);
  btas::Tie(ropr, btas::shape(0, 2), ropr_diag);

  btas::STArray<double,2> scr1;
  btas::Gemm(btas::NoTrans, btas::Trans, 1.0, lopr_diag, ropr_diag, 1.0, scr1);
  btas::Copy(scr1, diag, 1); // preserve quantum number of diag
}

} // namespace mpsxx

#endif // __MPSXX_COMPUTE_DIAGONAL_ELEMENTS_H
