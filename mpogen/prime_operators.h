/*!
 *  \file  prime_operator.h
 *  \brief Primitive site operator generator
 *
 *  \par Usage:
 *  prime_op_generator<mpsxx::fermionic::_cre_a>(op, l_index, r_index, 1.0);
 */

#ifndef _MPSXX_CXX11_PRIME_OPERATOR_H
#define _MPSXX_CXX11_PRIME_OPERATOR_H 1

#include <btas/TVector.h>
#include <btas/QSPARSE/QSDArray.h>

#include "Fermion/Quantum.h"

namespace mpsxx     {

namespace fermionic {

//! Physical index, i.e. |0>, |u>, |d>, |ud>
enum PHYSICAL_INDEX { vacuum = 0, alpha = 1, beta = 2, pair = 3 };

inline btas::Qshapes<Quantum> fock()
{
  btas::Qshapes<Quantum> q(4);
  q[vacuum] = Quantum(0, 0);
  q[alpha ] = Quantum(1, 1);
  q[beta  ] = Quantum(1,-1);
  q[pair  ] = Quantum(2, 0);
  return std::move(q);
}

//! Non-zero component of primitive operator
struct prime_op_component {
  //! Constructor
  prime_op_component
  (const PHYSICAL_INDEX& _bra, const PHYSICAL_INDEX& _ket, const double& _data) : bra(_bra), ket(_ket), data(_data) { }
  //! Return true if bra index has odd particle number
  bool parity() const { return (bra == alpha || bra == beta); }

  PHYSICAL_INDEX
    bra;   //!< bra index
  PHYSICAL_INDEX
    ket;   //!< ket index
  double
    data;  //!< non-zero element
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Static objects for fermionic operators
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! indentity
struct _identity {
  const static Quantum q() { return Quantum( 0,  0); }
  const static btas::TVector<prime_op_component, 4> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha)
struct _cre_a {
  const static Quantum q() { return Quantum(+1, +1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! destruct (alpha)
struct _des_a {
  const static Quantum q() { return Quantum(-1, -1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! create (beta)
struct _cre_b {
  const static Quantum q() { return Quantum(+1, -1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! destruct (beta)
struct _des_b {
  const static Quantum q() { return Quantum(-1, +1); }
  const static btas::TVector<prime_op_component, 2> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta)
struct _cre_a_cre_b {
  const static Quantum q() { return Quantum(+2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x create (alpha): this is taken to be normal order
struct _cre_b_cre_a {
  const static Quantum q() { return Quantum(+2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! destruct (alpha) x destruct (beta): this is taken to be normal order
struct _des_a_des_b {
  const static Quantum q() { return Quantum(-2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! destruct (beta) x destruct (alpha)
struct _des_b_des_a {
  const static Quantum q() { return Quantum(-2,  0); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (alpha) x destruct (alpha): counting operator (alpha)
struct _cre_a_des_a {
  const static Quantum q() { return Quantum( 0,  0); }
  const static btas::TVector<prime_op_component, 2> elements;
};

//! create (alpha) x destruct (beta): spin up operator
struct _cre_a_des_b {
  const static Quantum q() { return Quantum( 0, +2); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x destruct (alpha): spin down operator
struct _cre_b_des_a {
  const static Quantum q() { return Quantum( 0, -2); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x destruct (beta): counting operator (beta)
struct _cre_b_des_b {
  const static Quantum q() { return Quantum( 0,  0); }
  const static btas::TVector<prime_op_component, 2> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta) x destruct (beta)
struct _cre_a_cre_b_des_b {
  const static Quantum q() { return Quantum(+1, +1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x create (alpha) x destruct (alpha)
struct _cre_b_cre_a_des_a {
  const static Quantum q() { return Quantum(+1, -1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (alpha) x destruct (alpha) x destruct (beta)
struct _cre_a_des_a_des_b {
  const static Quantum q() { return Quantum(-1, +1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

//! create (beta) x destruct (beta) x destruct (alpha)
struct _cre_b_des_b_des_a {
  const static Quantum q() { return Quantum(-1, -1); }
  const static btas::TVector<prime_op_component, 1> elements;
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x destruct (alpha) x create (beta) x destruct (beta): referred to chemist's notation
struct _cre_a_des_a_cre_b_des_b {
  const static Quantum q() { return Quantum( 0,  0); }
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
(btas::QSDArray<4, Quantum>& op, const int& l_index, const int& r_index, const double& scale = 1.0)
{
  const Quantum& l_quanta = op.qshape(0)[l_index];
  const Quantum& r_quanta = op.qshape(3)[r_index];

//std::cout << "DEBUG[prime_op_generator]: q(L)::" << l_quanta << " + q(R)::" << r_quanta << " =? (-) q(S)::" << prime_op_obj::q() << std::endl;
  assert((r_quanta + l_quanta) == -prime_op_obj::q());
  if(std::fabs(scale) < 1.0e-20) return;

  btas::DArray<4> block(1,1,1,1);
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
(btas::QSDArray<4, mpsxx::Quantum>& op, const int& l_index, const int& r_index, const double& scale = 1.0)
{
  const Quantum& l_quanta = op.qshape(0)[l_index];
  const Quantum& r_quanta = op.qshape(3)[r_index];

  assert((r_quanta + l_quanta) == -prime_op_obj::q());
  if(std::fabs(scale) < 1.0e-20) return;

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

}; // fermionic

}; // mpsxx

#endif // _MPSXX_CXX11_PRIME_OPERATOR_H
