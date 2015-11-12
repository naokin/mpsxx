#include <iostream>
#include <iomanip>

#include "gen_qcmpo_utils.h"

// Right-Operator Table (swap_comp_block == false)
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

/// Generate product operators, O_ls = O_l x O_s, to construct non-zero MPO elements.
/// The dimension is exactly the same as DMRG complementary operator formalism.
std::vector<mpsxx::DMRG_OpType> mpsxx::gen_product_operators (
  const mpsxx::DMRG_OpType& l_op,
  const mpsxx::DMRG_OpType& s_op,
  const std::vector<int>& comp_index, bool swap_comp_block)
{
//std::cout << "DEBUG :: gen_product_operators :: 01 :: trying " << l_op.label() << " x " << s_op.label() << std::endl;
  std::vector<DMRG_OpType> ls_ops;

  op::QC::CATEGORY r_cat = op::QC::MultTable[l_op.category][s_op.category];
//std::cout << "DEBUG :: gen_product_operators :: 02 :: prefix " << op::QC::LabTable[r_cat] << std::endl;

  if(swap_comp_block) r_cat = op::QC::ConjTable[r_cat];

  int ls_index[4] = { -1, -1, -1, -1 };

  size_t n = 0;
  for(size_t i = 0; i < 2 && l_op.index[i] >= 0; ++i)
    ls_index[n++] = l_op.index[i];
  for(size_t i = 0; i < 2 && s_op.index[i] >= 0; ++i)
    ls_index[n++] = s_op.index[i];

//std::cout << "DEBUG :: gen_product_operators :: 03 :: n = " << n << std::endl;
  switch(r_cat) {
    // I(L) x Cij(S) x Ckl(R)
    // Ci(L) x Cj(S) x Ckl(R)
    // Cij(L) x I(S) x Ckl(R)
    case op::QC::CREA_CREB:
    case op::QC::CREB_CREA:
    case op::QC::CREA_DESA:
    case op::QC::CREA_DESB:
    case op::QC::CREB_DESA:
    case op::QC::CREB_DESB:
    case op::QC::DESA_CREA:
    case op::QC::DESA_CREB:
    case op::QC::DESB_CREA:
    case op::QC::DESB_CREB:
    case op::QC::DESA_DESB:
    case op::QC::DESB_DESA:
      // if Cii(R) is non-zero
      if(swap_comp_block) {
//std::cout << "DEBUG :: gen_product_operators :: 04 " << std::endl;
        for(size_t i = 0; i < comp_index.size(); ++i) {
          int ix = comp_index[i];
          ls_ops.push_back(DMRG_OpType(r_cat,ix,ix));
        }
      }
    case op::QC::CREA_CREA:
    case op::QC::CREB_CREB:
    case op::QC::DESA_DESA:
    case op::QC::DESB_DESB:
      if(swap_comp_block) {
//std::cout << "DEBUG :: gen_product_operators :: 05 " << std::endl;
        for(size_t i = 0; i < comp_index.size(); ++i) {
          int ix = comp_index[i];
          for(size_t j = 0; j < i; ++j) {
            int jx = comp_index[j];
            if(ix > jx)
              ls_ops.push_back(DMRG_OpType(r_cat,jx,ix));
            else
              ls_ops.push_back(DMRG_OpType(r_cat,ix,jx));
          }
        }
      }
      else {
//std::cout << "DEBUG :: gen_product_operators :: 06 " << std::endl;
        ls_ops.push_back(DMRG_OpType(r_cat,ls_index[0],ls_index[1]));
      }
      break;
    // I(L) x Ci(S) x Cjkl(R)
    // Ci(L) x I(S) x Cjkl(R)
    // I(L) x Cijk(S) x Cl(R)
    // Ci(L) x Cjk(S) x Cl(R)
    // Cij(L) x Ck(S) x Cl(R)
    // Cijk(L) x I(S) x Cl(R)
    case op::QC::CREA:
    case op::QC::CREB:
    case op::QC::DESA:
    case op::QC::DESB:
//std::cout << "DEBUG :: gen_product_operators :: 07 " << std::endl;
      ls_ops.push_back(DMRG_OpType(r_cat,ls_index[0]));
      break;
    case op::QC::CREA_COMP:
    case op::QC::CREB_COMP:
    case op::QC::DESA_COMP:
    case op::QC::DESB_COMP:
//std::cout << "DEBUG :: gen_product_operators :: 08 " << std::endl;
      for(size_t i = 0; i < comp_index.size(); ++i)
        ls_ops.push_back(DMRG_OpType(r_cat,comp_index[i]));
      break;
    // I(L) x I(S) x Cijkl(R)
    // I(L) x Cijkl(S) x I(R)
    // Ci(L) x Cjkl(S) x I(R)
    // Cij(L) x Ckl(S) x I(R)
    // Cijk(L) x Cl(S) x I(R)
    // Cijkl(L) x I(S) x I(R)
    case op::QC::H:
    case op::QC::I:
//std::cout << "DEBUG :: gen_product_operators :: 09 " << std::endl;
      break;
    default:
      abort();
  }
//std::cout << "DEBUG :: gen_product_operators :: 10 " << std::endl;

std::cout << "DEBUG :: gen_product_operators :: FF :: " << l_op << " x " << s_op << std::endl;
std::cout << "DEBUG :: gen_product_operators :: FF :: " << "--------------------------------------------------" << std::endl;
for(size_t i = 0; i < ls_ops.size(); ++i) {
std::cout << "DEBUG :: gen_product_operators :: FF :: " << ls_ops[i] << std::endl;
}
std::cout << std::endl;
  return ls_ops;
}
