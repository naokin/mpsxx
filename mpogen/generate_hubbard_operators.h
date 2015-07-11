#ifndef _MPSXX_CXX11_GENERATE_HUBBARD_OPERATORS_H
#define _MPSXX_CXX11_GENERATE_HUBBARD_OPERATORS_H 1

#include <string>

#include <legacy/DENSE/DArray.h>

#include <MPSblas.h>
#include <symmetry/Fermion/Quantum.h>

namespace mpsxx     {

namespace fermionic {

void generate_hubbard_operators
(MPO<double, Quantum>& mpos, const double& t = 1.0, const double& u = 1.0, const std::string& prefix = "./");

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_GENERATE_HUBBARD_OPERATORS_H
