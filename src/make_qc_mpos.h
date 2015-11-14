#ifndef __MPSXX_MAKE_QC_MPOS_H
#define __MPSXX_MAKE_QC_MPOS_H

#include <string>
#include <btas/DENSE/TArray.h>

namespace mpsxx {

int get_group (const size_t& Ndiv, const int& Norbs, const int& i, const int& j);

int get_group (const size_t& Ndiv, const int& Norbs, const int& i, const int& j, const int& k, const int& l);

void make_qc_mpos (
  const size_t& Ndivs,
  const size_t& Norbs,
  const double& Ecore,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint,
  const std::string& opname,
  const std::string& prefix,
  const bool& enable_swap_sweep,
  const bool& do_compress);

} // namespace mpsxx

#endif // __MPSXX_MAKE_QC_MPOS_H
