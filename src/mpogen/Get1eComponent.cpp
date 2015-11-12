#include "gen_qcmpo_utils.h"
#include <cassert>

/// Determine 1e-integral component of the product block: O(L) x O(S) x O(R)
double mpsxx::Get1eComponent (
  const mpsxx::DMRG_OpType& l_op,
  const mpsxx::DMRG_OpType& s_op,
  const mpsxx::DMRG_OpType& r_op,
  const OneIntArray& oneint)
{
  int idx1e[2] = { -1, -1 };

  op::RDM::CATEGORY cat1e;

  size_t n = 0;
  switch(l_op.category) {
    case op::QC::CREA:
    case op::QC::CREB:
    case op::QC::DESA:
    case op::QC::DESB:
      cat1e = (cat1e << (2u*n)) | op::QC2RDM[l_op.category];
      idx1e[n++] = l_op.index[0];
      break;
    default:
      break;
  }
  switch(s_op.category) {
    case op::QC::H:
      cat1e = (cat1e << (4u*n)) | (op::RDM::CREA << 2u) | op::RDM::DESA;
      idx1e[n++] = s_op.index[0];
      idx1e[n++] = s_op.index[0];
      break;
    case op::QC::CREA:
    case op::QC::CREB:
    case op::QC::DESA:
    case op::QC::DESB:
      cat1e = (cat1e << (2u*n)) | op::QC2RDM[s_op.category];
      idx1e[n++] = s_op.index[0];
      break;
    default:
      break;
  }
  switch(r_op.category) {
    case op::QC::CREA:
    case op::QC::CREB:
    case op::QC::DESA:
    case op::QC::DESB:
      cat1e = (cat1e << (2u*n)) | op::QC2RDM[r_op.category];
      idx1e[n++] = r_op.index[0];
      break;
    default:
      break;
  }
  assert(n == 2);

  double value = 0.0;
  switch(cat1e) {
    case 0x2: // 0010 = CREA_DESA
    case 0x7: // 0111 = CREB_DESB
      value = oneint(idx1e[0],idx1e[1]);
      break;
    case 0x8: // 1000 = DESA_CREA
    case 0xd: // 1101 = DESB_CREB
      value =-oneint(idx1e[0],idx1e[1]);
      break;
    default:
      break;
  }

  return value;
}
