
#include "integral_component.h"

/*! Determine 1-particle integral component for complementary operator construction in DMRG */
double mpsxx::fermionic::int1e_component
(const mpsxx::fermionic::BIT_OPERATOR_TYPE& l_op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& s_op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& r_op,
 const btas::DArray<2>& _int1e)
{
  using namespace bit_operator;
  btas::IVector<2> indxs;
  unsigned int itype = 0, n = 0;
  if((l_op & NORMAL & TYPE) == SINGLE) {
    indxs[n++] = (l_op & INDEX & FIRST) >>  INDEX_SHIFT;
    itype     |= (l_op & MASK  & FIRST) >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
  }
  if((s_op & NORMAL & TYPE) == SINGLE) {
    indxs[n++] = (s_op & INDEX & FIRST) >>  INDEX_SHIFT;
    itype     |= (s_op & MASK  & FIRST) >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
  }
  if((s_op & TYPE) == HAM) {
    // CreDes
    BIT_OPERATOR_TYPE op_index = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
    indxs[n++] = op_index;
    indxs[n++] = op_index;
    itype = 0x0000000d; // 11 01 = Cre(A)Des(A) [ = Cre(B)Des(B) ]
  }
  if((r_op & NORMAL & TYPE) == SINGLE) {
    indxs[n++] = (r_op & INDEX & FIRST) >>  INDEX_SHIFT;
    itype     |= (r_op & MASK  & FIRST) >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
  }
  assert(n == 2);

  size_t& ix = indxs[0];
  size_t& jx = indxs[1];

  double value = 0.0;
  // get commutation type
  if     (itype == 0x0000000d || itype == 0x00000008) // 1101 or 1000 : CreDes
    value =  _int1e(ix,jx);
  else if(itype == 0x00000007 || itype == 0x00000002) // 0111 or 0010 : DesCre
    value = -_int1e(ix,jx);
  return value;
}

/*! Determine 2-particle integral component for complementary operator construction in DMRG
 *
 *  \param l_op operator type for left  block
 *  \param s_op operator type for site  block
 *  \param r_op operator type for right block
 *  \param _int2e 2-particle integrals
 *
 *  Def. operator product
 *  : C[i,a] x C[j,a] -- 11 x 11 : F
 *  : C[i,a] x C[j,b] -- 11 x 10 : E
 *  : C[i,b] x C[j,a] -- 10 x 11 : B
 *  : C[i,b] x C[j,b] -- 10 x 10 : A
 *  : C[i,a] x D[j,a] -- 11 x 01 : D
 *  : C[i,a] x D[j,b] -- 11 x 00 : C
 *  : C[i,b] x D[j,a] -- 10 x 01 : 9
 *  : C[i,b] x D[j,b] -- 10 x 00 : 8
 *  : D[i,a] x C[j,a] -- 01 x 11 : 7
 *  : D[i,a] x C[j,b] -- 01 x 10 : 6
 *  : D[i,b] x C[j,a] -- 00 x 11 : 3
 *  : D[i,b] x C[j,b] -- 00 x 10 : 2
 *  : D[i,a] x D[j,a] -- 01 x 01 : 5
 *  : D[i,a] x D[j,b] -- 01 x 00 : 4
 *  : D[i,b] x D[j,a] -- 00 x 01 : 1
 *  : D[i,b] x D[j,b] -- 00 x 00 : 0
 *
 *  Independent integral contributions
 *
 * F: C[i,a] x C[j,a] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * D[k,a] x D[l,a] : 5
 * E: C[i,a] x C[j,b] -- ( - V[i,j,k,l] )              * D[k,a] x D[l,b] : 4
 *                    -- ( + V[i,j,l,k] )              * D[k,b] x D[l,a] : 1
 * B: C[i,b] x C[j,a] -- ( - V[i,j,k,l] )              * D[k,b] x D[l,a] : 1
 *                    -- ( + V[i,j,l,k] )              * D[k,a] x D[l,b] : 4
 * A: C[i,b] x C[j,b] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * D[k,b] x D[l,b] : 0
 *
 * D: C[i,a] x D[j,a] -- ( - V[i,j,l,k] + V[i,k,j,l] ) * C[k,a] x D[l,a] : D
 *                    -- ( + V[i,j,k,l] - V[i,k,j,l] ) * D[k,a] x C[l,a] : 7
 *                    -- ( + V[i,k,j,l] )              * C[k,b] x D[l,b] : 8
 *                    -- ( - V[i,k,j,l] )              * D[k,b] x C[l,b] : 2
 * C: C[i,a] x D[j,b] -- ( - V[i,j,l,k] )              * C[k,b] x D[l,a] : 9
 *                    -- ( + V[i,j,k,l] )              * D[k,a] x C[l,b] : 6
 * 9: C[i,b] x D[j,a] -- ( - V[i,j,l,k] )              * C[k,a] x D[l,b] : C
 *                    -- ( + V[i,j,k,l] )              * D[k,b] x C[l,a] : 3
 * 8: C[i,b] x D[j,b] -- ( - V[i,j,l,k] + V[i,k,j,l] ) * C[k,b] x D[l,b] : 8
 *                    -- ( + V[i,j,k,l] - V[i,k,j,l] ) * D[k,b] x C[l,b] : 2
 *                    -- ( + V[i,k,j,l] )              * C[k,a] x D[l,a] : D
 *                    -- ( - V[i,k,j,l] )              * D[k,a] x C[l,a] : 7
 *
 * 7: D[i,a] x C[j,a] -- ( + V[i,j,k,l] - V[i,k,j,l] ) * C[k,a] x D[l,a] : D
 *                    -- ( - V[i,j,l,k] + V[i,k,j,l] ) * D[k,a] x C[l,a] : 7
 *                    -- ( - V[i,k,j,l] )              * C[k,b] x D[l,b] : 8
 *                    -- ( + V[i,k,j,l] )              * D[k,b] x C[l,b] : 2
 * 6: D[i,a] x C[j,b] -- ( + V[i,j,k,l] )              * C[k,a] x D[l,b] : C
 *                    -- ( - V[i,j,l,k] )              * D[k,b] x C[l,a] : 3
 * 3: D[i,b] x C[j,a] -- ( + V[i,j,k,l] )              * C[k,b] x D[l,a] : 9
 *                    -- ( - V[i,j,l,k] )              * D[k,a] x C[l,b] : 6
 * 2: D[i,b] x C[j,b] -- ( + V[i,j,k,l] - V[i,k,j,l] ) * C[k,b] x D[l,b] : 8
 *                    -- ( - V[i,j,l,k] + V[i,k,j,l] ) * D[k,b] x C[l,b] : 2
 *                    -- ( - V[i,k,j,l] )              * C[k,a] x D[l,a] : D
 *                    -- ( + V[i,k,j,l] )              * D[k,a] x C[l,a] : 7
 *
 * 5: D[i,a] x D[j,a] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * C[k,a] x C[l,a] : F
 * 4: D[i,a] x D[j,b] -- ( - V[i,j,k,l] )              * C[k,a] x C[l,b] : E
 *                    -- ( + V[i,j,l,k] )              * C[k,b] x C[l,a] : B
 * 1: D[i,b] x D[j,a] -- ( - V[i,j,k,l] )              * C[k,b] x C[l,a] : B
 *                    -- ( + V[i,j,l,k] )              * C[k,a] x C[l,b] : E
 * 0: D[i,b] x D[j,b] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * C[k,b] x C[l,b] : A
 */
double mpsxx::fermionic::int2e_component
(const mpsxx::fermionic::BIT_OPERATOR_TYPE& l_op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& s_op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& r_op,
 const btas::DArray<4>& _int2e)
{
  using namespace bit_operator;
  /*! Table of integral type (256 = 4 x 4 x 4 x 4 tensor)
   * 
   *       1:   int2e(i,j,l,k)
   *       2:   int2e(i,j,k,l)
   *       3:   int2e(i,k,j,l)
   *       5: - int2e(i,j,l,k)
   *       6: - int2e(i,j,k,l)
   *       7: - int2e(i,k,j,l)
   *       9:   int2e(i,j,l,k) - int2e(i,j,k,l)
   *      10:   int2e(i,j,k,l) - int2e(i,k,j,l)
   *      11:   int2e(i,k,j,l) - int2e(i,j,l,k)
   */
  const static unsigned int m_int2e_table[256] = {
  /*---------------------------------------------------------------------*/
  /*   **  0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */
  /*---------------------------------------------------------------------*/
  /* 0 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9,  0,  0,  0,  0,  0,
  /* 1 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  1,  0,
  /* 2 */  0,  0, 11,  0,  0,  0,  0,  3, 10,  0,  0,  0,  0,  7,  0,  0,
  /* 3 */  0,  0,  0,  0,  0,  0,  5,  0,  0,  2,  0,  0,  0,  0,  0,  0,
  /* 4 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  6,  0,
  /* 5 */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9,
  /* 6 */  0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,
  /* 7 */  0,  0,  3,  0,  0,  0,  0, 11,  7,  0,  0,  0,  0, 10,  0,  0,
  /* 8 */  0,  0, 10,  0,  0,  0,  0,  7, 11,  0,  0,  0,  0,  3,  0,  0,
  /* 9 */  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  0,
  /* A */  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* B */  0,  6,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* C */  0,  0,  0,  0,  0,  0,  2,  0,  0,  5,  0,  0,  0,  0,  0,  0,
  /* D */  0,  0,  7,  0,  0,  0,  0, 10,  3,  0,  0,  0,  0, 11,  0,  0,
  /* E */  0,  1,  0,  0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* F */  0,  0,  0,  0,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
  };

  /*! Calculate integral type from operator types
   *
   *  C: [1*] / D: [0*]
   *  A: [*1] / B: [*0]
   */
  btas::IVector<4> indxs;
  unsigned int itype = 0, n = 0;
  BIT_OPERATOR_TYPE op_comp;
  BIT_OPERATOR_TYPE op_indx;
  switch(l_op & NORMAL & TYPE) {
    case SINGLE:
      indxs[n++] = (l_op & INDEX & FIRST) >> INDEX_SHIFT;
      itype |= (l_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      break;
    case DOUBLE:
      indxs[n++] = (l_op & INDEX & FIRST) >> INDEX_SHIFT;
      itype |= (l_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      indxs[n++] = (l_op & INDEX & SECOND);
      itype |= (l_op & MASK & SECOND) >>               (FIELD_SHIFT + 2*n - 8);
      break;
    default:
      break;
  }
  switch(s_op & TYPE) {
    case SINGLE:
      indxs[n++] = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      itype |= (s_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      break;
    case DOUBLE:
      indxs[n++] = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      itype |= (s_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      indxs[n++] = (s_op & INDEX & SECOND);
      itype |= (s_op & MASK & SECOND) >>               (FIELD_SHIFT + 2*n - 8);
      break;
    case SINGLE_COMP:
      op_comp = (s_op & MASK & FIRST) >> (INDEX_SHIFT + FIELD_SHIFT);
      op_indx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      // CreCreDes
      if(op_comp & 2) {
        indxs[n++] = op_indx; itype |=  (op_comp      << (8 - 2*n)); // Cre(S)
        indxs[n++] = op_indx; itype |= ((op_comp ^ 1) << (8 - 2*n)); // Cre(S*)
        indxs[n++] = op_indx; itype |= ((op_comp ^ 3) << (8 - 2*n)); // Des(S*)
      }
      // CreDesDes
      else {
        indxs[n++] = op_indx; itype |= ((op_comp ^ 3) << (8 - 2*n)); // Cre(S*)
        indxs[n++] = op_indx; itype |= ((op_comp ^ 1) << (8 - 2*n)); // Des(S*)
        indxs[n++] = op_indx; itype |=  (op_comp      << (8 - 2*n)); // Des(S)
      }
      break;
    case HAM:
      // CreDesCreDes
      op_indx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      indxs[n++] = op_indx;
      indxs[n++] = op_indx;
      indxs[n++] = op_indx;
      indxs[n++] = op_indx;
      itype = 0x000000d8; // 11 01 10 00 = Cre(A)Des(A)Cre(B)Des(B)
      break;
    default:
      break;
  }
  switch(r_op & NORMAL & TYPE) {
    case SINGLE:
      indxs[n++] = (r_op & INDEX & FIRST) >> INDEX_SHIFT;
      itype |= (r_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      break;
    case DOUBLE:
      indxs[n++] = (r_op & INDEX & FIRST) >> INDEX_SHIFT;
      itype |= (r_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      indxs[n++] = (r_op & INDEX & SECOND);
      itype |= (r_op & MASK & SECOND) >>               (FIELD_SHIFT + 2*n - 8);
      break;
    default:
      break;
  }
  assert(n == 4 && itype < 256);
  unsigned int m_type = m_int2e_table[itype];
//std::cout << "type = " << ((itype & 0x80) > 0) << ((itype & 0x40) > 0) << " "
//                       << ((itype & 0x20) > 0) << ((itype & 0x10) > 0) << " "
//                       << ((itype & 0x08) > 0) << ((itype & 0x04) > 0) << " "
//                       << ((itype & 0x02) > 0) << ((itype & 0x01) > 0) << " ["
//                       << std::setw(2) << m_type << "] :: " << std::flush;

  /*! Return scaling factor for complementary operator
   * 
   *  m_type:
   *       1:   int2e(i,j,l,k)
   *       2:   int2e(i,j,k,l)
   *       3:   int2e(i,k,j,l)
   *       5: - int2e(i,j,l,k)
   *       6: - int2e(i,j,k,l)
   *       7: - int2e(i,k,j,l)
   *       9:   int2e(i,j,l,k) - int2e(i,j,k,l)
   *      10:   int2e(i,j,k,l) - int2e(i,k,j,l)
   *      11:   int2e(i,k,j,l) - int2e(i,j,l,k)
   */
  size_t& ix = indxs[0];
  size_t& jx = indxs[1];
  size_t& kx = indxs[2];
  size_t& lx = indxs[3];
//std::cout << "[" << ix << "," << jx << "," << kx << "," << lx << "] :: " << std::flush;

  double value = 0.0;
  // get commutation type
  switch(m_type % 4) {
    case 1:
      value += _int2e(ix,jx,lx,kx); if(m_type & 8) value -= _int2e(ix,jx,kx,lx);
      break;
    case 2:
      value += _int2e(ix,jx,kx,lx); if(m_type & 8) value -= _int2e(ix,kx,jx,lx);
      break;
    case 3:
      value += _int2e(ix,kx,jx,lx); if(m_type & 8) value -= _int2e(ix,jx,lx,kx);
      break;
    default:
      break;
  }
  // get scaling +/-
  return (m_type & 4) ? -value : value;
}
