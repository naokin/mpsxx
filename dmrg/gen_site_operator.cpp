
#include "gen_site_operator.h"
#include "integral_component.h"
#include "prime_operators.h"

/*!
 *  Generate matrix product operator on single site
 */
void mpsxx::gen_site_operator (
        btas::QSTArray<double,4,fermion>& op,
  const mpsxx::mpogen::BIT_OPERATOR_TYPE& l_op, const size_t& l_index,
  const mpsxx::mpogen::BIT_OPERATOR_TYPE& s_op,
  const mpsxx::mpogen::BIT_OPERATOR_TYPE& r_op, const size_t& r_index,
  const double& Ecore,
  const btas::TArray<double,2>& oneint, const btas::TArray<double,4>& twoint)
{
  using namespace mpogen::bit_operator;
  double v_int1e = 0.0;
  double v_int2e = 0.0;
  if((s_op & TYPE) == IDEN) {
    // SWAP SWEEP:
    // CiCi x I -> DkDk * Vijkl
    // DiDi <- I x CkCk * Vijkl
    v_int2e = 1.0;
    if((l_op & TYPE) == DOUBLE && (r_op & TYPE) == DOUBLE && (l_op & INDEX) != (r_op & INDEX)) {
      v_int2e = mpogen::int2e_component(l_op, s_op, r_op, twoint);
    }
    mpogen::prime_op_generator<mpogen::_identity>(op, l_index, r_index, v_int2e);
  }
  else if((s_op & TYPE) == HAM) {
    // Ecore is added into H @ site 0
    if(!(s_op & INDEX) && fabs(Ecore) > 1.0e-16) {
      mpogen::prime_op_generator<mpogen::_identity>(op, l_index, r_index, Ecore);
    }
    // I x H -> H
    // H <- H x I
    v_int1e = mpogen::int1e_component(l_op, s_op, r_op, oneint);
    v_int2e = mpogen::int2e_component(l_op, s_op, r_op, twoint);
    mpogen::prime_op_generator<mpogen::_cre_a_des_a>            (op, l_index, r_index, v_int1e);
    mpogen::prime_op_generator<mpogen::_cre_b_des_b>            (op, l_index, r_index, v_int1e);
    mpogen::prime_op_generator<mpogen::_cre_a_des_a_cre_b_des_b>(op, l_index, r_index, v_int2e);
  }
  else if((s_op & TYPE) == SINGLE) {
    // FORWARD:
    //   Ci x Cj -> CiCj
    // CiCi x Dj -> Qk
    //
    // SWAP SWEEP:
    //   Ci x Cj -> DkDk
    // CiCi x Dj -> Dk
    v_int2e = 1.0;
    if     ((l_op & NORMAL & TYPE) == SINGLE && (r_op & TYPE) == DOUBLE) {
      mpogen::BIT_OPERATOR_TYPE conj = (l_op & COMP) ? (COMP | CONJ_S) : ZERO;
      if((DOUBLE | ((l_op ^ conj) & FIRST) | (s_op & FIRST) >> INDEX_SHIFT) != r_op
      && (DOUBLE | ((l_op ^ conj) & FIRST) >> INDEX_SHIFT | (s_op & FIRST)) != r_op) {
        v_int2e = mpogen::int2e_component((l_op ^ conj), s_op, r_op, twoint);
      }
    }
//#CURRENT LINE -- something wrong in this expression
    else if((l_op & TYPE) == DOUBLE && (r_op & NORMAL & TYPE) == SINGLE) {
      mpogen::BIT_OPERATOR_TYPE conj = (r_op & COMP) ? (COMP | CONJ_S) : ZERO;
      if((DOUBLE | (s_op & FIRST) | ((r_op ^ conj) & FIRST) >> INDEX_SHIFT) != l_op
      && (DOUBLE | (s_op & FIRST) >> INDEX_SHIFT | ((r_op ^ conj) & FIRST)) != l_op) {
        v_int2e = mpogen::int2e_component(l_op, s_op, (r_op ^ conj), twoint);
      }
    }
    switch(s_op & MASK) {
      case CRE_A:
        mpogen::prime_op_generator<mpogen::_cre_a>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B:
        mpogen::prime_op_generator<mpogen::_cre_b>(op, l_index, r_index, v_int2e);
        break;
      case DES_A:
        mpogen::prime_op_generator<mpogen::_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_B:
        mpogen::prime_op_generator<mpogen::_des_b>(op, l_index, r_index, v_int2e);
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
      mpogen::BIT_OPERATOR_TYPE conj = (l_op & COMP) ? (COMP | CONJ_S) : ZERO;
      v_int1e = mpogen::int1e_component((l_op ^ conj), s_op, r_op, oneint) * 0.5;
      v_int2e = mpogen::int2e_component((l_op ^ conj), s_op, r_op, twoint);
    }
    else if((r_op & NORMAL & TYPE) == SINGLE) {
      mpogen::BIT_OPERATOR_TYPE conj = (r_op & COMP) ? (COMP | CONJ_S) : ZERO;
      v_int1e = mpogen::int1e_component(l_op, s_op, (r_op ^ conj), oneint) * 0.5;
      v_int2e = mpogen::int2e_component(l_op, s_op, (r_op ^ conj), twoint);
    }
    switch(s_op & MASK) {
      case CRE_A:
        mpogen::prime_op_generator<mpogen::_cre_a>            (op, l_index, r_index, v_int1e);
        mpogen::prime_op_generator<mpogen::_cre_a_cre_b_des_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B:
        mpogen::prime_op_generator<mpogen::_cre_b>            (op, l_index, r_index, v_int1e);
        mpogen::prime_op_generator<mpogen::_cre_b_cre_a_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_A:
        mpogen::prime_op_generator<mpogen::_des_a>            (op, l_index, r_index, v_int1e);
        mpogen::prime_op_generator<mpogen::_cre_b_des_b_des_a>(op, l_index, r_index, v_int2e);
        break;
      case DES_B:
        mpogen::prime_op_generator<mpogen::_des_b>            (op, l_index, r_index, v_int1e);
        mpogen::prime_op_generator<mpogen::_cre_a_des_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        break;
    }
  }
  else if((s_op & TYPE) == DOUBLE) {
    v_int2e = 1.0;
    if     ((l_op & TYPE) == SINGLE && (r_op & TYPE) == SINGLE) { // swap sweep dir.
      v_int2e = mpogen::int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((l_op & TYPE) == SINGLE && (r_op & TYPE) == SINGLE_COMP) {
      v_int2e = mpogen::int2e_component(l_op, s_op, (r_op ^ COMP ^ CONJ_S), twoint);
    }
    else if((r_op & TYPE) == SINGLE && (l_op & TYPE) == SINGLE_COMP) {
      v_int2e = mpogen::int2e_component((l_op ^ COMP ^ CONJ_S), s_op, r_op, twoint);
    }
    else if((l_op & TYPE) == DOUBLE && (l_op & INDEX) != (s_op & INDEX)) {
      v_int2e = mpogen::int2e_component(l_op, s_op, r_op, twoint);
    }
    else if((r_op & TYPE) == DOUBLE && (r_op & INDEX) != (s_op & INDEX)) {
      v_int2e = mpogen::int2e_component(l_op, s_op, r_op, twoint);
    }
    switch(s_op & MASK) {
      case CRE_A_DES_A:
        mpogen::prime_op_generator<mpogen::_cre_a_des_a>(op, l_index, r_index, v_int2e);
        break;
      case CRE_B_DES_B:
        mpogen::prime_op_generator<mpogen::_cre_b_des_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_A_CRE_B:
        mpogen::prime_op_generator<mpogen::_cre_a_cre_b>(op, l_index, r_index, v_int2e);
        break;
      case CRE_A_DES_B:
        mpogen::prime_op_generator<mpogen::_cre_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      case DES_A_CRE_B:
//      mpogen::prime_op_generator<mpogen::_des_a_cre_b>(op, l_index, r_index, v_int2e);
        mpogen::prime_op_generator<mpogen::_cre_b_des_a>(op, l_index, r_index,-v_int2e);
        break;
      case DES_A_DES_B:
        mpogen::prime_op_generator<mpogen::_des_a_des_b>(op, l_index, r_index, v_int2e);
        break;
      default:
        // No other double operators on single site
        break;
    }
  }
//std::cout << "DEBUG :: " << mpogen::translate(l_op) << " x " << mpogen::translate(s_op) << " = " << mpogen::translate(r_op) << " with ";
//std::cout << "v1 = " << std::setprecision(2) << std::setw(5) << v_int1e << " :: "
//          << "v2 = " << std::setprecision(2) << std::setw(5) << v_int2e << " :: " << std::endl;
  return;
}

