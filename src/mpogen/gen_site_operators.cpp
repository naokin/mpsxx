#include "gen_qcmpo_utils.h"

/// Generate matrix product operator element on a single site
void mpsxx::gen_site_operators (
//      mpsxx::MPO<double,fermion>& mpo,
        btas::QSTArray<double,4,fermion>& mpo,
  const mpsxx::DMRG_OpType& l_op, size_t l_idx,
  const mpsxx::DMRG_OpType& s_op,
  const mpsxx::DMRG_OpType& r_op, size_t r_idx,
  const OneIntArray& oneint,
  const TwoIntArray& twoint)
{
  double val1e = 0.0;
  double val2e = 0.0;

  switch(s_op.category) {
    case op::QC::I:
      // Normal case
      // O(L) x I(S) -> O(R)
      val2e = 1.0;
      // Swap case
      // Cij(L) x I(S) -> Ckl(R) x Vijkl
      if(l_op.index[0] != r_op.index[0]
      && l_op.index[1] != r_op.index[1]) {
        val2e = Get2eComponent(l_op,s_op,r_op,twoint);
      }
      add_operators<op::QC::I>(op,l_idx,r_idx,val2e);
      break;
    case op::QC::H:
      // I x H -> H
      // H <- H x I
      val1e = Get1eComponent(l_op,s_op,r_op,oneint);
      val2e = Get2eComponent(l_op,s_op,r_op,twoint);
      add_operators<op::QC::CREA_DESA>(op,l_idx,r_idx,val1e);
      add_operators<op::QC::CREB_DESB>(op,l_idx,r_idx,val1e);
      add_operators<op::QC::CREA_DESA_CREB_DESB>(op,l_idx,r_idx,val2e);
      break;
    case op::QC::CREA_COMP:
    case op::QC::CREB_COMP:
    case op::QC::DESA_COMP:
    case op::QC::DESB_COMP:
      val2e = Get2eComponent(l_op,s_op,r_op,twoint);
  }
  else if((s_op & TYPE) == SINGLE) {
    // FORWARD:
    //   Ci x Cj -> CiCj
    // CiCi x Dj -> Qk
    //
    // SWAP SWEEP:
    //   Ci x Cj -> DkDk
    // CiCi x Dj -> Dk
    val2e = 1.0;
    if     ((l_op & NORMAL & TYPE) == SINGLE && (r_op & TYPE) == DOUBLE) {
      BIT_OPERATOR_TYPE conj = (l_op & COMP) ? (COMP | CONJ_S) : ZERO;
      if((DOUBLE | (l_op ^ conj) & FIRST | (s_op & FIRST) >> INDEX_SHIFT) != r_op) {
        val2e = int2e_component((l_op ^ conj), s_op, r_op, twoint);
      }
    }
    else if((l_op & TYPE) == DOUBLE && (r_op & NORMAL & TYPE) == SINGLE) {
      BIT_OPERATOR_TYPE conj = (r_op & COMP) ? (COMP | CONJ_S) : ZERO;
      if((DOUBLE | s_op & FIRST | ((r_op ^ conj) & FIRST) >> INDEX_SHIFT) != l_op) {
        val2e = int2e_component(l_op, s_op, (r_op ^ conj), twoint);
      }
    }
    switch(s_op & MASK) {
      case CRE_A:
        prime_op_generator<_cre_a>(op, l_index, r_index, val2e);
        break;
      case CRE_B:
        prime_op_generator<_cre_b>(op, l_index, r_index, val2e);
        break;
      case DES_A:
        prime_op_generator<_des_a>(op, l_index, r_index, val2e);
        break;
      case DES_B:
        prime_op_generator<_des_b>(op, l_index, r_index, val2e);
        break;
      default:
        break;
    }
  }
  else if((s_op & TYPE) == SINGLE_COMP) {
    // FORWARD:
    // Ci x Qi -> H
    // I  x Qk -> Qk
    //
    // SWAP SWEEP:
    // Ci x Qi -> I
    // I  x Qk -> Dk
    if     ((l_op & NORMAL & TYPE) == SINGLE) {
      BIT_OPERATOR_TYPE conj = (l_op & COMP) ? (COMP | CONJ_S) : ZERO;
      val1e = int1e_component((l_op ^ conj), s_op, r_op, oneint) * 0.5;
      val2e = int2e_component((l_op ^ conj), s_op, r_op, twoint);
    }
    else if((r_op & NORMAL & TYPE) == SINGLE) {
      BIT_OPERATOR_TYPE conj = (r_op & COMP) ? (COMP | CONJ_S) : ZERO;
      val1e = int1e_component(l_op, s_op, (r_op ^ conj), oneint) * 0.5;
      val2e = int2e_component(l_op, s_op, (r_op ^ conj), twoint);
    }
    switch(s_op & MASK) {
      case CRE_A:
        prime_op_generator<_cre_a>            (op, l_index, r_index, val1e);
        prime_op_generator<_cre_a_cre_b_des_b>(op, l_index, r_index, val2e);
        break;
      case CRE_B:
        prime_op_generator<_cre_b>            (op, l_index, r_index, val1e);
        prime_op_generator<_cre_b_cre_a_des_a>(op, l_index, r_index, val2e);
        break;
      case DES_A:
        prime_op_generator<_des_a>            (op, l_index, r_index, val1e);
        prime_op_generator<_cre_b_des_b_des_a>(op, l_index, r_index, val2e);
        break;
      case DES_B:
        prime_op_generator<_des_b>            (op, l_index, r_index, val1e);
        prime_op_generator<_cre_a_des_a_des_b>(op, l_index, r_index, val2e);
        break;
      default:
        break;
    }
  }
  else if((s_op & TYPE) == DOUBLE) {
    val2e = 1.0;
    if     ((l_op & TYPE) == SINGLE && (r_op & TYPE) == SINGLE) { // swap sweep dir.
      val2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((l_op & TYPE) == SINGLE && (r_op & TYPE) == SINGLE_COMP) {
      val2e = int2e_component(l_op, s_op, (r_op ^ COMP ^ CONJ_S), twoint);
    }
    else if((r_op & TYPE) == SINGLE && (l_op & TYPE) == SINGLE_COMP) {
      val2e = int2e_component((l_op ^ COMP ^ CONJ_S), s_op, r_op, twoint);
    }
    else if((l_op & TYPE) == DOUBLE && (l_op & INDEX) != (s_op & INDEX)) {
      val2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((r_op & TYPE) == DOUBLE && (r_op & INDEX) != (s_op & INDEX)) {
      val2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    switch(s_op & MASK) {
      case CRE_A_DES_A:
        prime_op_generator<_cre_a_des_a>(op, l_index, r_index, val2e);
        break;
      case CRE_B_DES_B:
        prime_op_generator<_cre_b_des_b>(op, l_index, r_index, val2e);
        break;
      case CRE_A_CRE_B:
        prime_op_generator<_cre_a_cre_b>(op, l_index, r_index, val2e);
        break;
      case CRE_A_DES_B:
        prime_op_generator<_cre_a_des_b>(op, l_index, r_index, val2e);
        break;
      case DES_A_CRE_B:
//      prime_op_generator<_des_a_cre_b>(op, l_index, r_index, val2e);
        prime_op_generator<_cre_b_des_a>(op, l_index, r_index,-val2e);
        break;
      case DES_A_DES_B:
        prime_op_generator<_des_a_des_b>(op, l_index, r_index, val2e);
        break;
      default:
        // No other double operators on single site
        break;
    }
  }
//std::cout << "v1 = " << std::setprecision(2) << std::setw(5) << val1e << " :: "
//          << "v2 = " << std::setprecision(2) << std::setw(5) << val2e << " :: " << std::flush;
  return;
}

