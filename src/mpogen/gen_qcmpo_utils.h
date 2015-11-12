#ifndef __MPSXX_MPOGEN_GEN_QCMPO_UTILS_H
#define __MPSXX_MPOGEN_GEN_QCMPO_UTILS_H

#include <vector>
#include <legacy/DENSE/TArray.h>

#include "DMRG_OpType.h"

namespace mpsxx {

typedef btas::TArray<double,2> OneIntArray;
typedef btas::TArray<double,4> TwoIntArray;

std::vector<DMRG_OpType> gen_product_operators (
  const DMRG_OpType& l_op,
  const DMRG_OpType& s_op,
  const std::vector<int>& comp_index, bool swap_comp_block = false);

//void gen_site_operators (
////      MPO<double,fermion>& mpo,
//        btas::QSTArray<double,4,fermion>& mpo,
//  const DMRG_OpType& l_op, size_t l_idx,
//  const DMRG_OpType& s_op,
//  const DMRG_OpType& r_op, size_t r_idx,
//  const OneIntArray& oneint,
//  const TwoIntArray& twoint);

double Get1eComponent (
  const DMRG_OpType& l_op,
  const DMRG_OpType& s_op,
  const DMRG_OpType& r_op,
  const OneIntArray& oneint);

double Get2eComponent (
  const DMRG_OpType& l_op,
  const DMRG_OpType& s_op,
  const DMRG_OpType& r_op,
  const TwoIntArray& twoint);

} // namespace mpsxx

#endif // __MPSXX_MPOGEN_GEN_QCMPO_UTILS_H
