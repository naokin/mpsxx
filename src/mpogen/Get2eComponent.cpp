#include "gen_qcmpo_utils.h"

/// Determine 2e-integral component of the product block: O(L) x O(S) x O(R)
/// \param l_op operator type for left  block
/// \param s_op operator type for site  block
/// \param r_op operator type for right block
/// \param twoint 2e integrals
///
///  Def. operator product
///  : C[i,a] x C[j,a] -- 11 x 11 : F
///  : C[i,a] x C[j,b] -- 11 x 10 : E
///  : C[i,b] x C[j,a] -- 10 x 11 : B
///  : C[i,b] x C[j,b] -- 10 x 10 : A
///  : C[i,a] x D[j,a] -- 11 x 01 : D
///  : C[i,a] x D[j,b] -- 11 x 00 : C
///  : C[i,b] x D[j,a] -- 10 x 01 : 9
///  : C[i,b] x D[j,b] -- 10 x 00 : 8
///  : D[i,a] x C[j,a] -- 01 x 11 : 7
///  : D[i,a] x C[j,b] -- 01 x 10 : 6
///  : D[i,b] x C[j,a] -- 00 x 11 : 3
///  : D[i,b] x C[j,b] -- 00 x 10 : 2
///  : D[i,a] x D[j,a] -- 01 x 01 : 5
///  : D[i,a] x D[j,b] -- 01 x 00 : 4
///  : D[i,b] x D[j,a] -- 00 x 01 : 1
///  : D[i,b] x D[j,b] -- 00 x 00 : 0
///
///  Independent integral contributions
///
/// F: C[i,a] x C[j,a] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * D[k,a] x D[l,a] : 5
/// E: C[i,a] x C[j,b] -- ( - V[i,j,k,l] )              * D[k,a] x D[l,b] : 4
///                    -- ( + V[i,j,l,k] )              * D[k,b] x D[l,a] : 1
/// B: C[i,b] x C[j,a] -- ( - V[i,j,k,l] )              * D[k,b] x D[l,a] : 1
///                    -- ( + V[i,j,l,k] )              * D[k,a] x D[l,b] : 4
/// A: C[i,b] x C[j,b] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * D[k,b] x D[l,b] : 0
///
/// D: C[i,a] x D[j,a] -- ( - V[i,j,l,k] + V[i,k,j,l] ) * C[k,a] x D[l,a] : D
///                    -- ( + V[i,j,k,l] - V[i,k,j,l] ) * D[k,a] x C[l,a] : 7
///                    -- ( + V[i,k,j,l] )              * C[k,b] x D[l,b] : 8
///                    -- ( - V[i,k,j,l] )              * D[k,b] x C[l,b] : 2
/// C: C[i,a] x D[j,b] -- ( - V[i,j,l,k] )              * C[k,b] x D[l,a] : 9
///                    -- ( + V[i,j,k,l] )              * D[k,a] x C[l,b] : 6
/// 9: C[i,b] x D[j,a] -- ( - V[i,j,l,k] )              * C[k,a] x D[l,b] : C
///                    -- ( + V[i,j,k,l] )              * D[k,b] x C[l,a] : 3
/// 8: C[i,b] x D[j,b] -- ( - V[i,j,l,k] + V[i,k,j,l] ) * C[k,b] x D[l,b] : 8
///                    -- ( + V[i,j,k,l] - V[i,k,j,l] ) * D[k,b] x C[l,b] : 2
///                    -- ( + V[i,k,j,l] )              * C[k,a] x D[l,a] : D
///                    -- ( - V[i,k,j,l] )              * D[k,a] x C[l,a] : 7
///
/// 7: D[i,a] x C[j,a] -- ( + V[i,j,k,l] - V[i,k,j,l] ) * C[k,a] x D[l,a] : D
///                    -- ( - V[i,j,l,k] + V[i,k,j,l] ) * D[k,a] x C[l,a] : 7
///                    -- ( - V[i,k,j,l] )              * C[k,b] x D[l,b] : 8
///                    -- ( + V[i,k,j,l] )              * D[k,b] x C[l,b] : 2
/// 6: D[i,a] x C[j,b] -- ( + V[i,j,k,l] )              * C[k,a] x D[l,b] : C
///                    -- ( - V[i,j,l,k] )              * D[k,b] x C[l,a] : 3
/// 3: D[i,b] x C[j,a] -- ( + V[i,j,k,l] )              * C[k,b] x D[l,a] : 9
///                    -- ( - V[i,j,l,k] )              * D[k,a] x C[l,b] : 6
/// 2: D[i,b] x C[j,b] -- ( + V[i,j,k,l] - V[i,k,j,l] ) * C[k,b] x D[l,b] : 8
///                    -- ( - V[i,j,l,k] + V[i,k,j,l] ) * D[k,b] x C[l,b] : 2
///                    -- ( - V[i,k,j,l] )              * C[k,a] x D[l,a] : D
///                    -- ( + V[i,k,j,l] )              * D[k,a] x C[l,a] : 7
///
/// 5: D[i,a] x D[j,a] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * C[k,a] x C[l,a] : F
/// 4: D[i,a] x D[j,b] -- ( - V[i,j,k,l] )              * C[k,a] x C[l,b] : E
///                    -- ( + V[i,j,l,k] )              * C[k,b] x C[l,a] : B
/// 1: D[i,b] x D[j,a] -- ( - V[i,j,k,l] )              * C[k,b] x C[l,a] : B
///                    -- ( + V[i,j,l,k] )              * C[k,a] x C[l,b] : E
/// 0: D[i,b] x D[j,b] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * C[k,b] x C[l,b] : A
double mpsxx::Get2eComponent (
  const mpsxx::DMRG_OpType& l_op,
  const mpsxx::DMRG_OpType& s_op,
  const mpsxx::DMRG_OpType& r_op,
  const TwoIntArray& twoint)
{
  /// Table of integral type (256 = 4 x 4 x 4 x 4 tensor)
  /// 
  ///      1:   int2e(i,j,l,k)
  ///      2:   int2e(i,j,k,l)
  ///      3:   int2e(i,k,j,l)
  ///      5: - int2e(i,j,l,k)
  ///      6: - int2e(i,j,k,l)
  ///      7: - int2e(i,k,j,l)
  ///      9:   int2e(i,j,l,k) - int2e(i,j,k,l)
  ///     10:   int2e(i,j,k,l) - int2e(i,k,j,l)
  ///     11:   int2e(i,k,j,l) - int2e(i,j,l,k)
  const static unsigned int Int2eTable[256] = {
  /*------------------------------------------------------------------------------------*/
  /*   **  0    1    2    3    4    5    6    7    8    9    A    B    C    D    E    F */
  /*------------------------------------------------------------------------------------*/
  /* 0 */  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  9u,  0u,  0u,  0u,  0u,  0u,
  /* 1 */  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  6u,  0u,  0u,  1u,  0u,
  /* 2 */  0u,  0u, 11u,  0u,  0u,  0u,  0u,  3u, 10u,  0u,  0u,  0u,  0u,  7u,  0u,  0u,
  /* 3 */  0u,  0u,  0u,  0u,  0u,  0u,  5u,  0u,  0u,  2u,  0u,  0u,  0u,  0u,  0u,  0u,
  /* 4 */  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  1u,  0u,  0u,  6u,  0u,
  /* 5 */  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  9u,
  /* 6 */  0u,  0u,  0u,  5u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  2u,  0u,  0u,  0u,
  /* 7 */  0u,  0u,  3u,  0u,  0u,  0u,  0u, 11u,  7u,  0u,  0u,  0u,  0u, 10u,  0u,  0u,
  /* 8 */  0u,  0u, 10u,  0u,  0u,  0u,  0u,  7u, 11u,  0u,  0u,  0u,  0u,  3u,  0u,  0u,
  /* 9 */  0u,  0u,  0u,  2u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  5u,  0u,  0u,  0u,
  /* A */  9u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,
  /* B */  0u,  6u,  0u,  0u,  1u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,
  /* C */  0u,  0u,  0u,  0u,  0u,  0u,  2u,  0u,  0u,  5u,  0u,  0u,  0u,  0u,  0u,  0u,
  /* D */  0u,  0u,  7u,  0u,  0u,  0u,  0u, 10u,  3u,  0u,  0u,  0u,  0u, 11u,  0u,  0u,
  /* E */  0u,  1u,  0u,  0u,  6u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,
  /* F */  0u,  0u,  0u,  0u,  0u,  9u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u,  0u
  };

  unsigned int idx2e[4];

  /// Calculate integral type from operator types
  /// CREA = [11] = 3
  /// CREB = [10] = 2
  /// DESA = [01] = 1
  /// DESB = [00] = 0
  unsigned int n = 0u;
  unsigned int l = 0u;
  switch(l_op.category) {
    case op::QC::CREA: ++l;
    case op::QC::CREB: ++l;
    case op::QC::DESA: ++l;
    case op::QC::DESB:
      idx2e[n++] = l_op.index[0];
      break;
    case op::QC::CREA_CREA: ++l;
    case op::QC::CREA_CREB: ++l;
    case op::QC::CREA_DESA: ++l;
    case op::QC::CREA_DESB: ++l;
    case op::QC::CREB_CREA: ++l;
    case op::QC::CREB_CREB: ++l;
    case op::QC::CREB_DESA: ++l;
    case op::QC::CREB_DESB: ++l;
    case op::QC::DESA_CREA: ++l;
    case op::QC::DESA_CREB: ++l;
    case op::QC::DESA_DESA: ++l;
    case op::QC::DESA_DESB: ++l;
    case op::QC::DESB_CREA: ++l;
    case op::QC::DESB_CREB: ++l;
    case op::QC::DESB_DESA: ++l;
    case op::QC::DESB_DESB:
      idx2e[n++] = l_op.index[0];
      idx2e[n++] = l_op.index[1];
      break;
    default:
      break;
  }
  unsigned int ls = 0u;
  switch(s_op.category) {
    case op::QC::H: // i = CREA_DESA_CREB_DESB = [1101 1000] = 16*13+ 8 = 216
      idx2e[n++] = s_op.index[0];
      l <<= 2u;
      ls +=160;
    case op::QC::CREA_COMP: // i = CREA_CREB_DESB = [11 1000] = 16*3+ 8 =  56
      ls +=  4;
    case op::QC::DESB_COMP: // i = CREA_DESA_DESB = [11 0100] = 16*3+ 4 =  52
      ls +=  7;
    case op::QC::CREB_COMP: // i = CREB_CREA_DESA = [10 1101] = 16*2+13 =  45
      ls += 12;
    case op::QC::DESA_COMP: // i = CREB_DESB_DESA = [10 0001] = 16*2+ 1 =  33
      idx2e[n++] = s_op.index[0];
      idx2e[n++] = s_op.index[0];
      idx2e[n++] = s_op.index[0];
      l <<= 6u;
      ls += 33;
      break;
    case op::QC::CREA_CREA: ++ls;
    case op::QC::CREA_CREB: ++ls;
    case op::QC::CREA_DESA: ++ls;
    case op::QC::CREA_DESB: ++ls;
    case op::QC::CREB_CREA: ++ls;
    case op::QC::CREB_CREB: ++ls;
    case op::QC::CREB_DESA: ++ls;
    case op::QC::CREB_DESB: ++ls;
    case op::QC::DESA_CREA: ++ls;
    case op::QC::DESA_CREB: ++ls;
    case op::QC::DESA_DESA: ++ls;
    case op::QC::DESA_DESB: ++ls;
    case op::QC::DESB_CREA: ++ls;
    case op::QC::DESB_CREB: ++ls;
    case op::QC::DESB_DESA: ++ls;
    case op::QC::DESB_DESB:
      idx2e[n++] = s_op.index[0];
      idx2e[n++] = s_op.index[1];
      l <<= 4u;
      break;
    default:
      break;
  }
  ls += l;
  unsigned int lsr = 0u;
  switch(r_op.category) {
    case op::QC::CREA: ++lsr;
    case op::QC::CREB: ++lsr;
    case op::QC::DESA: ++lsr;
    case op::QC::DESB:
      idx2e[n++] = r_op.index[0];
      ls <<= 2u;
      break;
    case op::QC::CREA_CREA: ++lsr;
    case op::QC::CREA_CREB: ++lsr;
    case op::QC::CREA_DESA: ++lsr;
    case op::QC::CREA_DESB: ++lsr;
    case op::QC::CREB_CREA: ++lsr;
    case op::QC::CREB_CREB: ++lsr;
    case op::QC::CREB_DESA: ++lsr;
    case op::QC::CREB_DESB: ++lsr;
    case op::QC::DESA_CREA: ++lsr;
    case op::QC::DESA_CREB: ++lsr;
    case op::QC::DESA_DESA: ++lsr;
    case op::QC::DESA_DESB: ++lsr;
    case op::QC::DESB_CREA: ++lsr;
    case op::QC::DESB_CREB: ++lsr;
    case op::QC::DESB_DESA: ++lsr;
    case op::QC::DESB_DESB:
      idx2e[n++] = r_op.index[0];
      idx2e[n++] = r_op.index[1];
      ls <<= 4u;
      break;
    default:
      break;
  }
  lsr += ls;

  assert(n == 4u && lsr < 256u);

  unsigned int t = Int2eTable[lsr];
  // DEBUG
  std::cout << "t = " << ((lsr & 0x80) > 0) << ((lsr & 0x40) > 0) << " "
                      << ((lsr & 0x20) > 0) << ((lsr & 0x10) > 0) << " "
                      << ((lsr & 0x08) > 0) << ((lsr & 0x04) > 0) << " "
                      << ((lsr & 0x02) > 0) << ((lsr & 0x01) > 0) << " ["
                      << std::setw(2) << t << "] :: ";

  /// Return scaling factor for complementary operator
  /// case t:
  ///      1:   int2e(i,j,l,k)
  ///      2:   int2e(i,j,k,l)
  ///      3:   int2e(i,k,j,l)
  ///      5: - int2e(i,j,l,k)
  ///      6: - int2e(i,j,k,l)
  ///      7: - int2e(i,k,j,l)
  ///      9:   int2e(i,j,l,k) - int2e(i,j,k,l)
  ///     10:   int2e(i,j,k,l) - int2e(i,k,j,l)
  ///     11:   int2e(i,k,j,l) - int2e(i,j,l,k)
  unsigned int& ix = idx2e[0];
  unsigned int& jx = idx2e[1];
  unsigned int& kx = idx2e[2];
  unsigned int& lx = idx2e[3];

  // DEBUG
  std::cout << "[" << ix << "," << jx << "," << kx << "," << lx << "] :: ";

  double value = 0.0;
  // get commutation type
  switch(t % 4u) {
    case 1u:
      value += twoint(ix,jx,lx,kx); if(t & 8u) value -= twoint(ix,jx,kx,lx);
      break;
    case 2u:
      value += twoint(ix,jx,kx,lx); if(t & 8u) value -= twoint(ix,kx,jx,lx);
      break;
    case 3u:
      value += twoint(ix,kx,jx,lx); if(t & 8u) value -= twoint(ix,jx,lx,kx);
      break;
    default:
      break;
  }
  // DEBUG
  std::cout.precision(8);
  if(t & 4u) std::cout << std::fixed << std::setw(12) << -value << std::endl;
  else       std::cout << std::fixed << std::setw(12) <<  value << std::endl;
  // get scaling +/-
  return (t & 4u) ? -value : value;
}
