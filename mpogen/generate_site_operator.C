
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
//std::cout << "DEBUG[generate_site_operator]: 01 -- called" << std::endl;
  if((s_op & TYPE) == IDEN) {
//std::cout << "DEBUG[generate_site_operator]: 02" << std::endl;
    if((l_op & TYPE) == DOUBLE && (r_op & TYPE) == DOUBLE) {
//std::cout << "DEBUG[generate_site_operator]: 03" << std::endl;
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
//std::cout << "DEBUG[generate_site_operator]: 04" << std::endl;
      prime_op_generator<_identity>(op, l_index, r_index, v_int2e);
    }
    else {
//std::cout << "DEBUG[generate_site_operator]: 05" << std::endl;
      prime_op_generator<_identity>(op, l_index, r_index, 1.0);
    }
  }
  else if((s_op & TYPE) == HAM) {
//std::cout << "DEBUG[generate_site_operator]: 06" << std::endl;
    v_int1e = int1e_component(l_op, s_op, r_op, oneint);
    v_int2e = int2e_component(l_op, s_op, r_op, twoint);
//std::cout << "DEBUG[generate_site_operator]: 07" << std::endl;
    prime_op_generator<_cre_a_des_a>            (op, l_index, r_index, v_int1e);
//std::cout << "DEBUG[generate_site_operator]: 08" << std::endl;
    prime_op_generator<_cre_b_des_b>            (op, l_index, r_index, v_int1e);
//std::cout << "DEBUG[generate_site_operator]: 09" << std::endl;
    prime_op_generator<_cre_a_des_a_cre_b_des_b>(op, l_index, r_index, v_int2e);
  }
  else if((s_op & TYPE) == SINGLE) {
//std::cout << "DEBUG[generate_site_operator]: 10 -- L::" << translate(l_op) << " x S::" << translate(s_op) << " = R::" << translate(r_op) << std::endl;
    v_int2e = 1.0;
    if     ((l_op & TYPE) == SINGLE && (r_op & TYPE) == DOUBLE) {
      if((DOUBLE | l_op & FIRST | (s_op & FIRST) >> INDEX_SHIFT) != r_op) {
//std::cout << "DEBUG[generate_site_operator]: 10 -- Comp / Forward" << std::endl;
        v_int2e = int2e_component(l_op, s_op, r_op, twoint);
      }
    }
    else if((l_op & TYPE) == DOUBLE && (r_op & TYPE) == SINGLE) {
      if((DOUBLE | s_op & FIRST | (r_op & FIRST) >> INDEX_SHIFT) != l_op) {
//std::cout << "DEBUG[generate_site_operator]: 10 -- Comp / Backward" << std::endl;
        v_int2e = int2e_component(l_op, s_op, r_op, twoint);
      }
    }
//std::cout << "DEBUG[generate_site_operator]: 11" << std::endl;
    switch(s_op & MASK) {
      case CRE_A:
//std::cout << "DEBUG[generate_site_operator]: 12" << std::endl;
        prime_op_generator<_cre_a>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B:
//std::cout << "DEBUG[generate_site_operator]: 13" << std::endl;
        prime_op_generator<_cre_b>(op, l_index, r_index, v_int2e);
        break;
      case DES_A:
//std::cout << "DEBUG[generate_site_operator]: 14" << std::endl;
        prime_op_generator<_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_B:
//std::cout << "DEBUG[generate_site_operator]: 15" << std::endl;
        prime_op_generator<_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        break;
    }
//std::cout << "DEBUG[generate_site_operator]: 16" << std::endl;
  }
  else if((s_op & TYPE) == SINGLE_COMP) {
//std::cout << "DEBUG[generate_site_operator]: 17" << std::endl;
    size_t sx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[generate_site_operator]: 18" << std::endl;
    v_int1e = int1e_component(l_op, s_op, r_op, oneint);
//std::cout << "DEBUG[generate_site_operator]: 19" << std::endl;
    v_int2e = int2e_component(l_op, s_op, r_op, twoint);
//std::cout << "DEBUG[generate_site_operator]: 20" << std::endl;
    switch(s_op & MASK) {
      case CRE_A:
        prime_op_generator<_cre_a>            (op, l_index, r_index, v_int1e);
//std::cout << "DEBUG[generate_site_operator]: 21" << std::endl;
        prime_op_generator<_cre_a_cre_b_des_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B:
        prime_op_generator<_cre_b>            (op, l_index, r_index, v_int1e);
//std::cout << "DEBUG[generate_site_operator]: 22" << std::endl;
        prime_op_generator<_cre_b_cre_a_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_A:
        prime_op_generator<_des_a>            (op, l_index, r_index, v_int1e);
//std::cout << "DEBUG[generate_site_operator]: 23" << std::endl;
        prime_op_generator<_cre_b_des_b_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_B:
        prime_op_generator<_des_b>            (op, l_index, r_index, v_int1e);
//std::cout << "DEBUG[generate_site_operator]: 24" << std::endl;
        prime_op_generator<_cre_a_des_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        break;
    }
//std::cout << "DEBUG[generate_site_operator]: 25" << std::endl;
  }
  else if((s_op & TYPE) == DOUBLE) {
//std::cout << "DEBUG[generate_site_operator]: 26 -- L::" << translate(l_op) << " x S::" << translate(s_op) << " = R::" << translate(r_op) << std::endl;
    v_int2e = 1.0;
    if((l_op & TYPE) == SINGLE && (r_op & TYPE) == SINGLE_COMP
    || (r_op & TYPE) == SINGLE && (l_op & TYPE) == SINGLE_COMP) {
//std::cout << "DEBUG[generate_site_operator]: 26 -- Comp / ?" << std::endl;
      v_int2e = int2e_component(l_op, s_op, r_op, twoint);
    }
//std::cout << "DEBUG[generate_site_operator]: 27" << std::endl;
    switch(s_op & MASK) {
      case CRE_A_DES_A:
        prime_op_generator<_cre_a_des_a>(op, l_index, r_index, v_int2e);
//std::cout << "DEBUG[generate_site_operator]: 28" << std::endl;
        break;
      case CRE_B_DES_B:
        prime_op_generator<_cre_b_des_b>(op, l_index, r_index, v_int2e);
//std::cout << "DEBUG[generate_site_operator]: 29" << std::endl;
        break;
      case CRE_A_CRE_B:
        prime_op_generator<_cre_a_cre_b>(op, l_index, r_index, v_int2e);
//std::cout << "DEBUG[generate_site_operator]: 30" << std::endl;
        break;
      case CRE_A_DES_B:
        prime_op_generator<_cre_a_des_b>(op, l_index, r_index, v_int2e);
//std::cout << "DEBUG[generate_site_operator]: 31" << std::endl;
        break;
      case DES_A_CRE_B:
//      prime_op_generator<_des_a_cre_b>(op, l_index, r_index, v_int2e);
        prime_op_generator<_cre_b_des_a>(op, l_index, r_index,-v_int2e);
//std::cout << "DEBUG[generate_site_operator]: 32" << std::endl;
        break;
      case DES_A_DES_B:
        prime_op_generator<_des_a_des_b>(op, l_index, r_index, v_int2e);
//std::cout << "DEBUG[generate_site_operator]: 33" << std::endl;
        break;
      default:
        // No other double operators on single site
        break;
    }
//std::cout << "DEBUG[generate_site_operator]: 34" << std::endl;
  }
  return;
}

