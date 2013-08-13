#ifndef _MPSXX_CXX11_GENERATE_QC_OPERATORS_H
#define _MPSXX_CXX11_GENERATE_QC_OPERATORS_H 1

#include <btas/DENSE/DArray.h>

#include "MpOperators.h"

namespace mpsxx     {

namespace fermionic {

std::string get_mpofile(const std::string& prefix, const int& index);

void generate_qc_operators(MpOperators& mpos, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint, bool enable_swap_sweep_dir = false);

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_GENERATE_QC_OPERATORS_H
