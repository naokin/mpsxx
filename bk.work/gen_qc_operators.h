#ifndef __MPSXX_GEN_QC_OPERATORS_H
#define __MPSXX_GEN_QC_OPERATORS_H

#include <string>

#include <legacy/DENSE/TArray.h>

namespace mpsxx {

void gen_qc_operators (
  const size_t& N,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint,
  const std::string& opname,
  const std::string& prefix = ".");

} // namespace mpsxx

#endif // __MPSXX_MPOGEN_GENERATE_QC_OPERATORS_H
