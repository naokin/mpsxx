#ifndef _MPSXX_CXX11_GENERATE_HUBBARD_OPERATORS_H
#define _MPSXX_CXX11_GENERATE_HUBBARD_OPERATORS_H 1

#include <string>

#include <btas/DENSE/DArray.h>

#include <MpOperators.h>
#include <symmetry/Fermion/Quantum.h>

namespace mpsxx     {

namespace fermionic {

void generate_hubbard_operators
(MpOperators<Quantum>& mpos, const double& t = 1.0, const double& u = 1.0, const std::string& prefix = "./");

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_GENERATE_HUBBARD_OPERATORS_H
