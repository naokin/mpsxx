
#include "prime_operators.h"

using mpsxx::mpogen::prime_op_component;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Static objects for fermionic operators
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! indentity
const btas::TVector<prime_op_component, 4>
mpsxx::mpogen::_identity::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::vacuum, mpsxx::MpSite<fermion>::vacuum, 1.0),
  prime_op_component(mpsxx::MpSite<fermion>::alpha,  mpsxx::MpSite<fermion>::alpha,  1.0),
  prime_op_component(mpsxx::MpSite<fermion>::beta,   mpsxx::MpSite<fermion>::beta,   1.0),
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::pair,   1.0)
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha)
const btas::TVector<prime_op_component, 2>
mpsxx::mpogen::_cre_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::alpha,  mpsxx::MpSite<fermion>::vacuum, 1.0),
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::beta,  -1.0)
};

//! destruct (alpha)
const btas::TVector<prime_op_component, 2>
mpsxx::mpogen::_des_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::vacuum, mpsxx::MpSite<fermion>::alpha,  1.0),
  prime_op_component(mpsxx::MpSite<fermion>::beta,   mpsxx::MpSite<fermion>::pair,  -1.0)
};

//! create (beta)
const btas::TVector<prime_op_component, 2>
mpsxx::mpogen::_cre_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::beta,   mpsxx::MpSite<fermion>::vacuum, 1.0),
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::alpha,  1.0)
};

//! destruct (beta)
const btas::TVector<prime_op_component, 2>
mpsxx::mpogen::_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::vacuum, mpsxx::MpSite<fermion>::beta,   1.0),
  prime_op_component(mpsxx::MpSite<fermion>::alpha,  mpsxx::MpSite<fermion>::pair,   1.0)
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta)
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_a_cre_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::vacuum,-1.0)
};

//! create (beta) x create (alpha): this is taken to be normal order
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_b_cre_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::vacuum, 1.0)
};

//! destruct (alpha) x destruct (beta): this is taken to be normal order
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_des_a_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::vacuum, mpsxx::MpSite<fermion>::pair,   1.0)
};

//! destruct (beta) x destruct (alpha)
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_des_b_des_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::vacuum, mpsxx::MpSite<fermion>::pair,  -1.0)
};

//! create (alpha) x destruct (alpha): counting operator (alpha)
const btas::TVector<prime_op_component, 2>
mpsxx::mpogen::_cre_a_des_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::alpha,  mpsxx::MpSite<fermion>::alpha,  1.0),
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::pair,   1.0)
};

//! create (alpha) x destruct (beta): spin up operator
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_a_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::alpha,  mpsxx::MpSite<fermion>::beta,   1.0)
};

//! create (beta) x destruct (alpha): spin down operator
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_b_des_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::beta,   mpsxx::MpSite<fermion>::alpha,  1.0)
};

//! create (beta) x destruct (beta): counting operator (beta)
const btas::TVector<prime_op_component, 2>
mpsxx::mpogen::_cre_b_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::beta,   mpsxx::MpSite<fermion>::beta,   1.0),
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::pair,   1.0)
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta) x destruct (beta)
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_a_cre_b_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::beta,  -1.0)
};

//! create (beta) x create (alpha) x destruct (alpha)
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_b_cre_a_des_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::alpha,  1.0)
};

//! create (alpha) x destruct (alpha) x destruct (beta)
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_a_des_a_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::alpha,  mpsxx::MpSite<fermion>::pair,   1.0)
};

//! create (beta) x destruct (beta) x destruct (alpha)
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_b_des_b_des_a::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::beta,   mpsxx::MpSite<fermion>::pair,  -1.0)
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x destruct (alpha) x create (beta) x destruct (beta): referred to chemist's notation
const btas::TVector<prime_op_component, 1>
mpsxx::mpogen::_cre_a_des_a_cre_b_des_b::elements = {
  prime_op_component(mpsxx::MpSite<fermion>::pair,   mpsxx::MpSite<fermion>::pair,   1.0)
};

