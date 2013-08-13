
#include "get_product_ops.h"

// Site-Operator Table (sweep swap == false)
// -----+------------------------------+---------------------+
// L \ R| I    H    Ci   Di   Rl   Sl  | CiCj CiDj DiCj DiDj |
// -----+------------------------------+---------------------+
// I    | I    H    Cj   Dj   Ri   Si  | CjCj CjDj DjCj DjDj |
// H    | 0    I    0    0    0    0   | 0    0    0    0    |
// Ci   | 0    Dj   I    0    CjDj DjDj| Cj   Dj   0    0    |
// Di   | 0    Cj   0    I    CjCj CjDj| 0    0    Cj   Dj   |
// Rl   | 0    Dl   0    0    I    0   | 0    0    0    0    |
// Sl   | 0    Cl   0    0    0    I   | 0    0    0    0    |
// -----+------------------------------+---------------------+
// CiCj | 0    DkDk 0    0    Dk   0   | I    0    0    0    |
// CiDj | 0    CkDk 0    0    Ck   Dk  | 0    I    0    0    |
// DiCj | 0    CkDk 0    0    Ck   Dk  | 0    0    I    0    |
// DiDj | 0    CkCk 0    0    0    Ck  | 0    0    0    I    |
// -----+------------------------------+---------------------+

// Right-Operator Table (sweep swap == false)
// -----+------------------------------+---------------------+
// L \ S| I    H    Ck   Dk   Rl   Sl  | CkCk CkDk DkCk DkDk |
// -----+------------------------------+---------------------+
// I    | I    H    Ck   Dk   Rl   Sl  | CkCk CkDk DkCk DkDk |
// H    | H    0    0    0    0    0   | 0    0    0    0    |
// Ci   | Ci   0    CiCk CiDk 0    H   | 0    Rl   Rl   Sl   |
// Di   | Di   0    DiCk DiDk H    0   | Rl   Sl   Sl   0    |
// Rl   | Rl   0    0    H    0    0   | 0    0    0    0    |
// Sl   | Sl   0    H    0    0    0   | 0    0    0    0    |
// -----+------------------------------+---------------------+
// CiCj | CiCj 0    0    Rl   0    0   | 0    0    0    H    |
// CiDj | CiDj 0    Rl   Sl   0    0   | 0    H    H    0    |
// DiCj | DiCj 0    Rl   Sl   0    0   | 0    H    H    0    |
// DiDj | DiDj 0    Sl   0    0    0   | H    0    0    0    |
// -----+------------------------------+---------------------+

//! Generate product operator: O_ls = O_l x O_s
/*!
 *  Make non-zero MPO contributions of which the dimension is exactly the same as DMRG complementary operator formula
 *  MEMO: operator hierarchy: (high) [!, ~], [<<, >>], [==, !=], [&], [^], [|], [&&], [||] (low)
 *  TODO: This algorithm can be also useful in conventional DMRG code?
 */
std::vector<mpsxx::fermionic::BIT_OPERATOR_TYPE> mpsxx::fermionic::get_product_ops
(const BIT_OPERATOR_TYPE& l_op, const BIT_OPERATOR_TYPE& s_op, const std::vector<size_t>& r_indxs, bool _swap_sweep_dir)
{
  using namespace bit_operator;
  std::vector<BIT_OPERATOR_TYPE> ls_ops;
  // S == Identity
  if(s_op == IDEN) {
    // L = CiCj and R = DkDl, etc...
    if((l_op & NORMAL & TYPE) == DOUBLE) {
      if(_swap_sweep_dir) {
        // Fix operator type for right block
        BIT_OPERATOR_TYPE r_op = (l_op ^ CONJ_D) & MASK;
        // Check quantum number to have Pii or Qii component
        bool has_diag = ((r_op & FIRST) >> INDEX_SHIFT) != (r_op & SECOND);
        // Create product block operators
        r_op |= DOUBLE;
        for(size_t i = 0; i < r_indxs.size(); ++i) {
          size_t ix = r_indxs[i];
          if(has_diag)
            ls_ops.push_back(r_op | (ix << INDEX_SHIFT) | ix);
          for(size_t j = 0; j < i; ++j) {
            size_t jx = r_indxs[j];
            if(ix > jx)
              ls_ops.push_back(r_op | (jx << INDEX_SHIFT) | ix);
            else
              ls_ops.push_back(r_op | (ix << INDEX_SHIFT) | jx);
          }
        }
      }
      else {
        ls_ops.push_back(l_op);
      }
    }
    // L = Ci or CiComp -> R = DiComp or Di, etc...
    else if((l_op & NORMAL & TYPE) == SINGLE_1 && _swap_sweep_dir) {
      ls_ops.push_back(l_op ^ CONJ_S_1 ^ COMP);
    }
    else if(l_op & IDEN) {
      // Needs conjugation due to swap sweep direction?
      BIT_OPERATOR_TYPE conj = _swap_sweep_dir ? COMP : ZERO;
      ls_ops.push_back(l_op ^ conj);
    }
  }
  // L == Identity
  else if(l_op == IDEN) {
    // S = CiCi and R = DkDl, etc...
    if((s_op & NORMAL & TYPE) == DOUBLE) {
      if(_swap_sweep_dir) {
        // Fix operator type for right block
        BIT_OPERATOR_TYPE r_op = (s_op ^ CONJ_D) & MASK;
        // Check quantum number to have Pii or Qii component
        bool has_diag = ((r_op & FIRST) >> INDEX_SHIFT) != (r_op & SECOND);
        // Create product block operators
        for(size_t i = 0; i < r_indxs.size(); ++i) {
          size_t ix = r_indxs[i];
          if(has_diag)
            ls_ops.push_back(r_op | (ix << INDEX_SHIFT) | ix);
          for(size_t j = 0; j < i; ++j) {
            size_t jx = r_indxs[j];
            if(ix > jx)
              ls_ops.push_back(r_op | (jx << INDEX_SHIFT) | ix);
            else
              ls_ops.push_back(r_op | (ix << INDEX_SHIFT) | jx);
          }
        }
      }
      else {
        ls_ops.push_back(s_op);
      }
    }
    // S = CiCiDj -> Rk or S = Ci
    else if((s_op & NORMAL & TYPE) == SINGLE_1) {
      // Needs conjugation due to swap sweep direction?
      BIT_OPERATOR_TYPE conj = _swap_sweep_dir ? (CONJ_S_1 | COMP) : ZERO;
      // Create product block operators
      if(s_op & COMP) {
        // Fix operator type for right block
        BIT_OPERATOR_TYPE r_op = (s_op ^ conj) & (~INDEX);
        for(size_t i = 0; i < r_indxs.size(); ++i) {
          size_t ix = r_indxs[i];
          ls_ops.push_back(r_op | (ix << INDEX_SHIFT));
        }
      }
      else {
        ls_ops.push_back(s_op ^ conj);
      }
    }
    else if(s_op & IDEN) {
      // Needs conjugation due to swap sweep direction?
      BIT_OPERATOR_TYPE conj = _swap_sweep_dir ? COMP : ZERO;
      ls_ops.push_back((s_op & TYPE) ^ conj);
    }
  }
  // L = Ci or Rk and S = Cj or Rj
  else if((l_op & NORMAL & TYPE) == SINGLE_1 && (s_op & NORMAL & TYPE) == SINGLE_1) {
    // Ci(L) x Cj(S) -> CiCj(LS)
    if((l_op & COMP) == 0 && (s_op & COMP) == 0) {
      BIT_OPERATOR_TYPE r_op = DOUBLE | (l_op & FIRST) | ((s_op & FIRST) >> INDEX_SHIFT);
      if(_swap_sweep_dir) {
        // Fix operator type for right block
        r_op = (r_op ^ CONJ_D) & MASK;
        // Check quantum number to have Pii or Qii component
        bool has_diag = ((r_op & FIRST) >> INDEX_SHIFT) != (r_op & SECOND);
        // Create product block operators
        for(size_t i = 0; i < r_indxs.size(); ++i) {
          size_t ix = r_indxs[i];
          if(has_diag)
            ls_ops.push_back(r_op | (ix << INDEX_SHIFT) | ix);
          for(size_t j = 0; j < i; ++j) {
            size_t jx = r_indxs[j];
            if(ix > jx)
              ls_ops.push_back(r_op | (jx << INDEX_SHIFT) | ix);
            else
              ls_ops.push_back(r_op | (ix << INDEX_SHIFT) | jx);
          }
        }
      }
      else {
        ls_ops.push_back(r_op);
      }
    }
    // Rj(L) x Cj(S) -> H
    else if((l_op & COMP) && ((l_op ^ CONJ_S_1) & NORMAL) == s_op) {
      if(_swap_sweep_dir)
        ls_ops.push_back(IDEN);
      else
        ls_ops.push_back(HAM);
    }
    // Ci(L) x Rj(S) -> H
    else if((l_op & COMP) == 0 && ((s_op ^ CONJ_S_1) & MASK) == (l_op & MASK)) {
      if(_swap_sweep_dir)
        ls_ops.push_back(IDEN);
      else
        ls_ops.push_back(HAM);
    }
  }
  // L = CiCj and S = Dk, etc...
  else if((l_op & NORMAL & TYPE) == DOUBLE && (s_op & NORMAL & TYPE) == SINGLE_1 && (s_op & COMP) == 0) {
    // Checking quantum number
    Quantum q_ls = get_quantum(l_op) + get_quantum(s_op);
    if(abs(q_ls.p()) == 1 && abs(q_ls.s()) == 1) {
      // Needs conjugation due to swap sweep direction?
      BIT_OPERATOR_TYPE conj = _swap_sweep_dir ? (CONJ_S_1 | COMP) : ZERO;
      // Fix operator type for right block
      BIT_OPERATOR_TYPE r_op = (l_op & FIRST) ^ ((l_op & SECOND) << INDEX_SHIFT) ^ (s_op & FIRST);
                        r_op = (COMP | SINGLE_1 | (r_op & MASK & FIRST)) ^ conj;
      // Create product block operators
      for(size_t i = 0; i < r_indxs.size(); ++i) {
        size_t ix = r_indxs[i];
        ls_ops.push_back(r_op | (ix << INDEX_SHIFT));
      }
    }
  }
  // L = Ci and S = CjDk, etc...
  else if((s_op & NORMAL & TYPE) == DOUBLE && (l_op & NORMAL & TYPE) == SINGLE_1 && (l_op & COMP) == 0) {
    // Checking quantum number
    Quantum q_ls = get_quantum(l_op) + get_quantum(s_op);
    if(abs(q_ls.p()) == 1 && abs(q_ls.s()) == 1) {
      // Needs conjugation due to swap sweep direction?
      BIT_OPERATOR_TYPE conj = _swap_sweep_dir ? (CONJ_S_1 | COMP) : ZERO;
      // Fix operator type for right block
      BIT_OPERATOR_TYPE r_op = (s_op & FIRST) ^ ((s_op & SECOND) << INDEX_SHIFT) ^ (l_op & FIRST);
                        r_op = (COMP | SINGLE_1 | (r_op & MASK & FIRST)) ^ conj;
      // Create product block operators
      for(size_t i = 0; i < r_indxs.size(); ++i) {
        size_t ix = r_indxs[i];
        ls_ops.push_back(r_op | (ix << INDEX_SHIFT));
      }
    }
  }
  // L = CiCj and S = DkDl, etc...
  else if((l_op & NORMAL & TYPE) == DOUBLE && (s_op & NORMAL & TYPE) == DOUBLE) {
    // Checking quantum number
    Quantum q_ls = get_quantum(l_op) + get_quantum(s_op);
    if(q_ls == Quantum::zero()) {
      if(_swap_sweep_dir)
        ls_ops.push_back(IDEN);
      else
        ls_ops.push_back(HAM);
    }
  }

  return std::move(ls_ops);
}
