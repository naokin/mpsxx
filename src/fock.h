#ifndef _MPSXX_CXX11_FOCK_H
#define _MPSXX_CXX11_FOCK_H 1

#include <btas/DENSE/DArray.h>

namespace mpsxx {

void twopdm_to_onepdm
(const size_t& Norbs,
 const size_t& Nelec,
 const btas::DArray<4>& twopdm,
       btas::DArray<2>& onepdm);

void compute_fock_matrix
(const size_t& Norbs,
 const btas::DArray<2>& oneint,
 const btas::DArray<4>& twoint,
 const btas::DArray<2>& onepdm,
 const btas::DArray<4>& twopdm,
       btas::DArray<2>& fock);

void compute_gradient
(const size_t& Norbs,
 const btas::DArray<2>& fock,
       btas::DArray<2>& grad);

void compute_hessian
(const size_t& Norbs,
 const btas::DArray<2>& fock,
 const btas::DArray<2>& oneint,
 const btas::DArray<4>& twoint,
 const btas::DArray<2>& onepdm,
 const btas::DArray<4>& twopdm,
       btas::DArray<4>& hess);

double compute_energy
(const size_t& Norbs,
 const btas::DArray<2>& oneint,
 const btas::DArray<4>& twoint,
 const btas::DArray<2>& onepdm,
 const btas::DArray<4>& twopdm);

}; // namespace mpsxx

#endif // _MPSXX_CXX11_FOCK_H
