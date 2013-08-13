
#include "prime_operators.h"

using mpsxx::fermionic::prime_op_component;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Static objects for fermionic operators
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! indentity
const btas::TVector<prime_op_component, 4>
mpsxx::fermionic::_identity::elements = { prime_op_component(vacuum, vacuum, 1.0), prime_op_component(alpha, alpha, 1.0),
                                          prime_op_component(beta, beta, 1.0), prime_op_component(pair, pair, 1.0) };

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha)
const btas::TVector<prime_op_component, 2>
mpsxx::fermionic::_cre_a::elements = { prime_op_component(alpha, vacuum, 1.0), prime_op_component(pair, beta,-1.0) };

//! destruct (alpha)
const btas::TVector<prime_op_component, 2>
mpsxx::fermionic::_des_a::elements = { prime_op_component(vacuum, alpha, 1.0), prime_op_component(beta, pair,-1.0) };

//! create (beta)
const btas::TVector<prime_op_component, 2>
mpsxx::fermionic::_cre_b::elements = { prime_op_component(beta, vacuum, 1.0), prime_op_component(pair, alpha, 1.0) };

//! destruct (beta)
const btas::TVector<prime_op_component, 2>
mpsxx::fermionic::_des_b::elements = { prime_op_component(vacuum, beta, 1.0), prime_op_component(alpha, pair, 1.0) };

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta)
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_a_cre_b::elements = { prime_op_component(pair, vacuum,-1.0) };

//! create (beta) x create (alpha): this is taken to be normal order
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_b_cre_a::elements = { prime_op_component(pair, vacuum, 1.0) };

//! destruct (alpha) x destruct (beta): this is taken to be normal order
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_des_a_des_b::elements = { prime_op_component(vacuum, pair, 1.0) };

//! destruct (beta) x destruct (alpha)
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_des_b_des_a::elements = { prime_op_component(vacuum, pair,-1.0) };

//! create (alpha) x destruct (alpha): counting operator (alpha)
const btas::TVector<prime_op_component, 2>
mpsxx::fermionic::_cre_a_des_a::elements = { prime_op_component(alpha, alpha, 1.0), prime_op_component(pair, pair, 1.0) };

//! create (alpha) x destruct (beta): spin up operator
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_a_des_b::elements = { prime_op_component(alpha, beta, 1.0) };

//! create (beta) x destruct (alpha): spin down operator
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_b_des_a::elements = { prime_op_component(beta, alpha, 1.0) };

//! create (beta) x destruct (beta): counting operator (beta)
const btas::TVector<prime_op_component, 2>
mpsxx::fermionic::_cre_b_des_b::elements = { prime_op_component(beta, beta, 1.0), prime_op_component(pair, pair, 1.0) };

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x create (beta) x destruct (beta)
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_a_cre_b_des_b::elements = { prime_op_component(pair, beta,-1.0) };

//! create (beta) x create (alpha) x destruct (alpha)
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_b_cre_a_des_a::elements = { prime_op_component(pair, alpha, 1.0) };

//! create (alpha) x destruct (alpha) x destruct (beta)
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_a_des_a_des_b::elements = { prime_op_component(alpha, pair, 1.0) };

//! create (beta) x destruct (beta) x destruct (alpha)
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_b_des_b_des_a::elements = { prime_op_component(beta, pair,-1.0) };

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! create (alpha) x destruct (alpha) x create (beta) x destruct (beta): referred to chemist's notation
const btas::TVector<prime_op_component, 1>
mpsxx::fermionic::_cre_a_des_a_cre_b_des_b::elements = { prime_op_component(pair, pair, 1.0) };

