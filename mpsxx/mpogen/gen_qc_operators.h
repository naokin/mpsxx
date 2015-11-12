#ifndef __MPSXX_GEN_QC_OPERATORS_H
#define __MPSXX_GEN_QC_OPERATORS_H

#include <string>

#include <legacy/DENSE/TArray.h>

#include <MPX_Types.h>
#include <symmetry/fermion.h>

namespace mpsxx {

void gen_qc_operators (MPOs<double,fermion>& mpos, const btas::TArray<double,2>& oneint, const btas::TArray<double,4>& twoint, bool enable_swap_sweep_dir = false, const std::string& prefix = "./");

} // namespace mpsxx

#endif // ___MPSXX_MPOGEN_GENERATE_QC_OPERATORS_H
