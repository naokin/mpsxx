#ifndef _MPSXX_CXX11_DIAGONAL_ELEMENTS_H
#define _MPSXX_CXX11_DIAGONAL_ELEMENTS_H 1

#include <btas/SPARSE/SDdiagonal.h>
#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx {

//! Compute diagonal H elements for one-site algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSDArray<4, Q>& mpo0,
 const btas::QSDArray<3, Q>& lopr,
 const btas::QSDArray<3, Q>& ropr,
       btas::QSDArray<3, Q>& diag)
{
  btas::SDArray<3> mpo0_diag;
  btas::SDArray<2> lopr_diag;
  btas::SDArray<2> ropr_diag;

  btas::SDdiagonal(mpo0, btas::shape(1, 2), mpo0_diag);
  btas::SDdiagonal(lopr, btas::shape(0, 2), lopr_diag);
  btas::SDdiagonal(ropr, btas::shape(0, 2), ropr_diag);

  btas::SDArray<3> scr1;
  btas::SDcontract(1.0, lopr_diag, btas::shape(1), mpo0_diag, btas::shape(0), 1.0, scr1);
  btas::SDArray<3> scr2;
  btas::SDcontract(1.0, scr1,      btas::shape(2), ropr_diag, btas::shape(1), 1.0, scr2);
  btas::SDcopy(scr2, diag, 1); // preserve quantum number of diag
}

//! Compute diagonal H elements for two-site algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSDArray<4, Q>& lmpo,
 const btas::QSDArray<4, Q>& rmpo,
 const btas::QSDArray<3, Q>& lopr,
 const btas::QSDArray<3, Q>& ropr,
       btas::QSDArray<4, Q>& diag)
{
  btas::SDArray<3> lmpo_diag;
  btas::SDArray<3> rmpo_diag;
  btas::SDArray<2> lopr_diag;
  btas::SDArray<2> ropr_diag;

  btas::SDdiagonal(lmpo, btas::shape(1, 2), lmpo_diag);
  btas::SDdiagonal(rmpo, btas::shape(1, 2), rmpo_diag);
  btas::SDdiagonal(lopr, btas::shape(0, 2), lopr_diag);
  btas::SDdiagonal(ropr, btas::shape(0, 2), ropr_diag);

  btas::SDArray<3> scr1;
  btas::SDcontract(1.0, lopr_diag, btas::shape(1), lmpo_diag, btas::shape(0), 1.0, scr1);
  btas::SDArray<4> scr2;
  btas::SDcontract(1.0, scr1,      btas::shape(2), rmpo_diag, btas::shape(0), 1.0, scr2);
  btas::SDArray<4> scr3;
  btas::SDcontract(1.0, scr2,      btas::shape(3), ropr_diag, btas::shape(1), 1.0, scr3);
  btas::SDcopy(scr3, diag, 1); // preserve quantum number of diag
}

//! Compute diagonal H elements for merged algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSDArray<3, Q>& lopr,
 const btas::QSDArray<3, Q>& ropr,
       btas::QSDArray<2, Q>& diag)
{
  btas::SDArray<2> lopr_diag;
  btas::SDArray<2> ropr_diag;

  btas::SDdiagonal(lopr, btas::shape(0, 2), lopr_diag);
  btas::SDdiagonal(ropr, btas::shape(0, 2), ropr_diag);

  btas::SDArray<2> scr1;
  btas::SDgemm(btas::NoTrans, btas::Trans, 1.0, lopr_diag, ropr_diag, 1.0, scr1);
  btas::SDcopy(scr1, diag, 1); // preserve quantum number of diag
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_DIAGONAL_ELEMENTS_H
