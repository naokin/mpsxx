#ifndef __MPSXX_MPOGEN_GENERATE_HUBBARD_OPERATORS_H
#define __MPSXX_MPOGEN_GENERATE_HUBBARD_OPERATORS_H 1

#include <string>

#include <legacy/DENSE/TArray.h>

#include <MPX_types.h>
#include <symmetry/fermion.h>

namespace mpsxx     {

void gen_hubbard_operators (MPOs<double,fermion>& mpos, const double& t = 1.0, const double& u = 1.0, const std::string& prefix = "./");

}; // namespace mpsxx

#endif // __MPSXX_MPOGEN_GENERATE_HUBBARD_OPERATORS_H
