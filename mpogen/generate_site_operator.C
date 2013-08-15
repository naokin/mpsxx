
#include "generate_site_operator.h"
#include "integral_component.h"
#include "prime_operators.h"

/*!
 *  Generate matrix product operator on single site
 */
void mpsxx::fermionic::generate_site_operator
(btas::QSDArray<4, mpsxx::fermionic::Quantum>& op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& l_op, const size_t& l_index,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& s_op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& r_op, const size_t& r_index,
 const btas::DArray<2>& oneint, const btas::DArray<4>& twoint)
{
  using namespace bit_operator;
  double v_int1e = 0.0;
  double v_int2e = 0.0;
  if((s_op & TYPE) == IDEN) {
    v_int2e = 1.0;
    if((l_op & TYPE) == DOUBLE && (r_op & TYPE) == DOUBLE && (l_op & INDEX) != (r_op & INDEX)) {
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    prime_op_generator<_identity>(op, l_index, r_index, v_int2e);
  }
  else if((s_op & TYPE) == HAM) {
    v_int1e = int1e_component(l_op, s_op, r_op, oneint) * 0.5;
    v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    prime_op_generator<_cre_a_des_a>            (op, l_index, r_index, v_int1e);
    prime_op_generator<_cre_b_des_b>            (op, l_index, r_index, v_int1e);
    prime_op_generator<_cre_a_des_a_cre_b_des_b>(op, l_index, r_index, v_int2e);
  }
  else if((s_op & TYPE) == SINGLE) {
    v_int2e = 1.0;
//  if     ((l_op & TYPE) == SINGLE && (r_op & TYPE) == DOUBLE) {
    if     ((l_op & NORMAL & TYPE) == SINGLE && (r_op & TYPE) == DOUBLE) {
      BIT_OPERATOR_TYPE conj = (l_op & COMP) ? (COMP | CONJ_S) : ZERO;
      if((DOUBLE | (l_op ^ conj) & FIRST | (s_op & FIRST) >> INDEX_SHIFT) != r_op) {
        v_int2e = int2e_component((l_op ^ conj), s_op, r_op, twoint);
      }
    }
//  else if((l_op & TYPE) == DOUBLE && (r_op & TYPE) == SINGLE) {
    else if((l_op & TYPE) == DOUBLE && (r_op & NORMAL & TYPE) == SINGLE) {
      BIT_OPERATOR_TYPE conj = (r_op & COMP) ? (COMP | CONJ_S) : ZERO;
      if((DOUBLE | s_op & FIRST | ((r_op ^ conj) & FIRST) >> INDEX_SHIFT) != l_op) {
        v_int2e = int2e_component(l_op, s_op, (r_op ^ conj), twoint);
      }
    }
    switch(s_op & MASK) {
      case CRE_A:
        prime_op_generator<_cre_a>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B:
        prime_op_generator<_cre_b>(op, l_index, r_index, v_int2e);
        break;
      case DES_A:
        prime_op_generator<_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_B:
        prime_op_generator<_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        break;
    }
  }
  else if((s_op & TYPE) == SINGLE_COMP) {
    v_int1e = int1e_component(l_op, s_op, r_op, oneint) * 0.5;
    v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    switch(s_op & MASK) {
      case CRE_A:
        prime_op_generator<_cre_a>            (op, l_index, r_index, v_int1e);
        prime_op_generator<_cre_a_cre_b_des_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B:
        prime_op_generator<_cre_b>            (op, l_index, r_index, v_int1e);
        prime_op_generator<_cre_b_cre_a_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_A:
        prime_op_generator<_des_a>            (op, l_index, r_index, v_int1e);
        prime_op_generator<_cre_b_des_b_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_B:
        prime_op_generator<_des_b>            (op, l_index, r_index, v_int1e);
        prime_op_generator<_cre_a_des_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        break;
    }
  }
  else if((s_op & TYPE) == DOUBLE) {
    v_int2e = 1.0;
    if     ((l_op & TYPE) == SINGLE && (r_op & TYPE) == SINGLE_COMP) {
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((r_op & TYPE) == SINGLE && (l_op & TYPE) == SINGLE_COMP) {
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((l_op & TYPE) == DOUBLE && (l_op & INDEX) != (s_op & INDEX)) {
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((r_op & TYPE) == DOUBLE && (r_op & INDEX) != (s_op & INDEX)) {
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    }
    switch(s_op & MASK) {
      case CRE_A_DES_A:
        prime_op_generator<_cre_a_des_a>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B_DES_B:
        prime_op_generator<_cre_b_des_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_A_CRE_B:
        prime_op_generator<_cre_a_cre_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_A_DES_B:
        prime_op_generator<_cre_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      case DES_A_CRE_B:
//      prime_op_generator<_des_a_cre_b>(op, l_index, r_index, v_int2e);
        prime_op_generator<_cre_b_des_a>(op, l_index, r_index,-v_int2e);
        break;
      case DES_A_DES_B:
        prime_op_generator<_des_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        // No other double operators on single site
        break;
    }
  }
//std::cout << "v1 = " << v_int1e << ", v2 = " << v_int2e << " :: " << std::flush;
  return;
}

