/*!
 *  \file  prime_operator.h
 *  \brief Primitive site operator generator
 *
 *  \par Usage:
 *  prime_op_generator<mpsxx::mpogen::_cre_a>(op, l_index, r_index, 1.0);
 */

#ifndef __MPSXX_MPOGEN_PRIME_OPERATOR_H
#define __MPSXX_MPOGEN_PRIME_OPERATOR_H 1

//#include <btas/TVector.h>
#include <btas/QSPARSE/QSTArray.h>

#include "mpsite.h"
#include "fermion.h"

namespace mpsxx {
namespace mpogen {

//! Non-zero component of primitive operator
struct prime_op_component {
  //! Constructor
  prime_op_component (const MpSite<fermion>::PHYSICAL_INDEX& _bra, const MpSite<fermion>::PHYSICAL_INDEX& _ket, const double& _data) : bra(_bra), ket(_ket), data(_data) { }

  //! Return true if bra index has odd particle number
  bool parity() const { return (bra == MpSite<fermion>::alpha || bra == MpSite<fermion>::beta); }

  MpSite<fermion>::PHYSICAL_INDEX bra;   //!< bra index

  MpSite<fermion>::PHYSICAL_INDEX ket;   //!< ket index

  double data;  //!< non-zero element
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Static objects for fermionic operators
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! indentity
struct _identity {
  const static fermion q() { return fermion( 0,  0); }
  const static btas::TVector<prime_op_component, 4> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha)
struct _cre_a {
  const static fermion q() { return fermion(+1, +1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! destruct (alpha)
struct _des_a {
  const static fermion q() { return fermion(-1, -1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! create (beta)
struct _cre_b {
  const static fermion q() { return fermion(+1, -1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! destruct (beta)
struct _des_b {
  const static fermion q() { return fermion(-1, +1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta)
struct _cre_a_cre_b {
  const static fermion q() { return fermion(+2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x create (alpha): this is taken to be normal order
struct _cre_b_cre_a {
  const static fermion q() { return fermion(+2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! destruct (alpha) x destruct (beta): this is taken to be normal order
struct _des_a_des_b {
  const static fermion q() { return fermion(-2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! destruct (beta) x destruct (alpha)
struct _des_b_des_a {
  const static fermion q() { return fermion(-2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (alpha) x destruct (alpha): counting operator (alpha)
struct _cre_a_des_a {
  const static fermion q() { return fermion( 0,  0); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! create (alpha) x destruct (beta): spin up operator
struct _cre_a_des_b {
  const static fermion q() { return fermion( 0, +2); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x destruct (alpha): spin down operator
struct _cre_b_des_a {
  const static fermion q() { return fermion( 0, -2); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x destruct (beta): counting operator (beta)
struct _cre_b_des_b {
  const static fermion q() { return fermion( 0,  0); }
  const static btas::TVector<prime_op_component, 2> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta) x destruct (beta)
struct _cre_a_cre_b_des_b {
  const static fermion q() { return fermion(+1, +1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x create (alpha) x destruct (alpha)
struct _cre_b_cre_a_des_a {
  const static fermion q() { return fermion(+1, -1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (alpha) x destruct (alpha) x destruct (beta)
struct _cre_a_des_a_des_b {
  const static fermion q() { return fermion(-1, +1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x destruct (beta) x destruct (alpha)
struct _cre_b_des_b_des_a {
  const static fermion q() { return fermion(-1, -1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x destruct (alpha) x create (beta) x destruct (beta): referred to chemist's notation
struct _cre_a_des_a_cre_b_des_b {
  const static fermion q() { return fermion( 0,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! Prime operator generator
/*! generates primitive site mpo
 *
 *  \param prime_op_obj site operator type object, e.g. _cre_a, _cre_a_des_a, etc...
 *  \param l_index  left  index
 *  \param r_index  right index
 *  \param scale scaling factor, set 1. by default
 */
template<class prime_op_obj>
void prime_op_generator
(btas::QSTArray<double,4, fermion>& op, const int& l_index, const int& r_index, const double& scale = 1.0)
{
  const fermion& l_quanta = op.qshape(0)[l_index];
  const fermion& r_quanta = op.qshape(3)[r_index];

//std::cout << "DEBUG[prime_op_generator]: q(L)::" << l_quanta << " + q(R)::" << r_quanta << " =? (-) q(S)::" << prime_op_obj::q() << std::endl;
  assert((r_quanta * l_quanta) == -prime_op_obj::q());
  if(std::fabs(scale) < 1.0e-16) return;

  btas::TArray<double,4> block(1,1,1,1);
  for(const prime_op_component& p : prime_op_obj::elements) {
    if(l_quanta.parity() && p.parity())
      block =-scale * p.data;
    else
      block = scale * p.data;
    op.insert(btas::shape(l_index, p.bra, p.ket, r_index), block);
  }
}

//
// Real-sparse array version?
/*

template<class prime_op_obj>
void prime_op_generator
(btas::QSTArray<double,4, mpsxx::fermion>& op, const int& l_index, const int& r_index, const double& scale = 1.0)
{
  const fermion& l_quanta = op.qshape(0)[l_index];
  const fermion& r_quanta = op.qshape(3)[r_index];

  assert((r_quanta + l_quanta) == -prime_op_obj::q());
  if(std::fabs(scale) < 1.0e-16) return;

  double block;
  for(const prime_op_component& p : prime_op_obj::elements) {
    if(l_quanta.parity() && p.parity())
      block =-scale * p.data;
    else
      block = scale * p.data;
    op.insert(btas::shape(l_index, p.bra, p.ket, l_index), block);
  }
}

 */

} // mpogen
} // mpsxx

#endif // __MPSXX_MPOGEN_PRIME_OPERATOR_H
