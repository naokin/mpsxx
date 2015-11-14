#ifndef __MPSXX_MPOGEN_GENERATE_SITE_OPERATOR_H
#define __MPSXX_MPOGEN_GENERATE_SITE_OPERATOR_H 1

#include <btas/DENSE/TArray.h>
#include <btas/QSPARSE/QSTArray.h>

#include "bit_operator_type.h"

namespace mpsxx     {

/// Generate matrix product operator on single site
void gen_site_operator (
        btas::QSTArray<double,4,fermion>& op,
  const mpogen::BIT_OPERATOR_TYPE& l_op, const size_t& l_index,
  const mpogen::BIT_OPERATOR_TYPE& s_op,
  const mpogen::BIT_OPERATOR_TYPE& r_op, const size_t& r_index,
  const double& Ecore,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint);

} // namespace mpsxx

#endif // __MPSXX_MPOGEN_GENERATE_SITE_OP_H
