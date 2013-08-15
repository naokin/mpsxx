
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
  BIT_OPERATOR_TYPE op_indx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
  switch(l_op & NORMAL & TYPE) {
    case SINGLE:
      indxs[n++] = (l_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      itype     |= (l_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      break;
    case DOUBLE:
      indxs[n++] = (l_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      itype     |= (l_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      indxs[n++] = (l_op & INDEX & SECOND);
      itype     |= (l_op & MASK  & SECOND)         >>               (FIELD_SHIFT + 2*n - 4);
      break;
    default:
      break;
  }
  switch(s_op & TYPE) {
    case SINGLE:
      indxs[n++] = (s_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      itype     |= (s_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      break;
    case SINGLE_COMP:
      indxs[n++] = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      if((l_op & TYPE) == SINGLE_COMP)
        itype   |= (s_op & MASK  & FIRST ^ CONJ_S) >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      else
        itype   |= (s_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      break;
    case DOUBLE:
      indxs[n++] = (s_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      itype     |= (s_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      indxs[n++] = (s_op & INDEX & SECOND);
      itype     |= (s_op & MASK  & SECOND)         >>               (FIELD_SHIFT + 2*n - 4);
      break;
    case HAM:
      // CreDesCreDes
      op_indx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      for(; n < 2; n++) indxs[n] = op_indx;
      itype = 0x0000000d; // 11 01 = Cre(A)Des(A) [ = Cre(B)Des(B) ]
    default:
      break;
  }
  switch(r_op & TYPE) {
    case SINGLE:
      indxs[n++] = (r_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      itype     |= (r_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      break;
    case SINGLE_COMP:
      indxs[n++] = (r_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      if((s_op & TYPE) == SINGLE_COMP)
        itype   |= (r_op & MASK & FIRST ^ CONJ_S)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      else
        itype   |= (r_op & MASK & FIRST)           >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      break;
    case DOUBLE:
      indxs[n++] = (r_op & INDEX & FIRST)          >>  INDEX_SHIFT;
      itype     |= (r_op & MASK  & FIRST)          >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 4);
      indxs[n++] = (r_op & INDEX & SECOND);
      itype     |= (r_op & MASK  & SECOND)         >>               (FIELD_SHIFT + 2*n - 4);
      break;
    default:
      break;
  }
  assert(n == 2);

  int& ix = indxs[0];
  int& jx = indxs[1];

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
//std::cout << "DEBUG[int2e_component]: 01" << std::endl;
      indxs[n++] = (l_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 02" << std::endl;
      itype |= (l_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      break;
    case DOUBLE:
//std::cout << "DEBUG[int2e_component]: 03" << std::endl;
      indxs[n++] = (l_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 04" << std::endl;
      itype |= (l_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
//std::cout << "DEBUG[int2e_component]: 05" << std::endl;
      indxs[n++] = (l_op & INDEX & SECOND);
//std::cout << "DEBUG[int2e_component]: 06" << std::endl;
      itype |= (l_op & MASK & SECOND) >>               (FIELD_SHIFT + 2*n - 8);
      break;
    default:
      break;
  }
  switch(s_op & TYPE) {
    case SINGLE:
//std::cout << "DEBUG[int2e_component]: 07" << std::endl;
      indxs[n++] = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 08" << std::endl;
      itype |= (s_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      break;
    case DOUBLE:
//std::cout << "DEBUG[int2e_component]: 09" << std::endl;
      indxs[n++] = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 10 -- n::" << n << std::endl;
      itype |= (s_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
//std::cout << "DEBUG[int2e_component]: 11" << std::endl;
      indxs[n++] = (s_op & INDEX & SECOND);
//std::cout << "DEBUG[int2e_component]: 12 -- n::" << n << std::endl;
      itype |= (s_op & MASK & SECOND) >>               (FIELD_SHIFT + 2*n - 8);
      break;
    case SINGLE_COMP:
//std::cout << "DEBUG[int2e_component]: 13" << std::endl;
      op_comp = (s_op & MASK & FIRST) >> (INDEX_SHIFT + FIELD_SHIFT);
//std::cout << "DEBUG[int2e_component]: 14" << std::endl;
      op_indx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
      // CreCreDes
      if(op_comp & 2) {
//std::cout << "DEBUG[int2e_component]: 15" << std::endl;
        indxs[n++] = op_indx; itype |=  (op_comp      << (8 - 2*n)); // Cre(S)
//std::cout << "DEBUG[int2e_component]: 16" << std::endl;
        indxs[n++] = op_indx; itype |= ((op_comp ^ 1) << (8 - 2*n)); // Cre(S*)
//std::cout << "DEBUG[int2e_component]: 17" << std::endl;
        indxs[n++] = op_indx; itype |= ((op_comp ^ 3) << (8 - 2*n)); // Des(S*)
      }
      // CreDesDes
      else {
//std::cout << "DEBUG[int2e_component]: 18" << std::endl;
        indxs[n++] = op_indx; itype |= ((op_comp ^ 3) << (8 - 2*n)); // Cre(S*)
//std::cout << "DEBUG[int2e_component]: 19" << std::endl;
        indxs[n++] = op_indx; itype |= ((op_comp ^ 1) << (8 - 2*n)); // Des(S*)
//std::cout << "DEBUG[int2e_component]: 20" << std::endl;
        indxs[n++] = op_indx; itype |=  (op_comp      << (8 - 2*n)); // Des(S)
      }
      break;
    case HAM:
      // CreDesCreDes
//std::cout << "DEBUG[int2e_component]: 21" << std::endl;
      op_indx = (s_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 22" << std::endl;
      for(; n < 4; n++) indxs[n] = op_indx;
//std::cout << "DEBUG[int2e_component]: 23" << std::endl;
      itype = 0x000000d8; // 11 01 10 00 = Cre(A)Des(A)Cre(B)Des(B)
      break;
    default:
      break;
  }
  switch(r_op & NORMAL & TYPE) {
    case SINGLE:
//std::cout << "DEBUG[int2e_component]: 24" << std::endl;
      indxs[n++] = (r_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 25 -- n::" << n << std::endl;
      itype |= (r_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
      break;
    case DOUBLE:
//std::cout << "DEBUG[int2e_component]: 26" << std::endl;
      indxs[n++] = (r_op & INDEX & FIRST) >> INDEX_SHIFT;
//std::cout << "DEBUG[int2e_component]: 27 -- n::" << n << std::endl;
      itype |= (r_op & MASK & FIRST)  >> (INDEX_SHIFT + FIELD_SHIFT + 2*n - 8);
//std::cout << "DEBUG[int2e_component]: 28" << std::endl;
      indxs[n++] = (r_op & INDEX & SECOND);
//std::cout << "DEBUG[int2e_component]: 29 -- n::" << n << std::endl;
      itype |= (r_op & MASK & SECOND) >>               (FIELD_SHIFT + 2*n - 8);
      break;
    default:
      break;
  }
//std::cout << "DEBUG[int2e_component]: 30" << std::endl;
  assert(n == 4 && itype < 256);
  unsigned int m_type = m_int2e_table[itype];

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
  int& ix = indxs[0];
  int& jx = indxs[1];
  int& kx = indxs[2];
  int& lx = indxs[3];

  double value = 0.0;
  // get commutation type
  switch(m_type % 4) {
    case 1:
//std::cout << "DEBUG[int2e_component]: 31" << std::endl;
      value += _int2e(ix,jx,lx,kx); if(m_type & 8) value -= _int2e(ix,jx,kx,lx);
      break;
    case 2:
//std::cout << "DEBUG[int2e_component]: 32 -- ix::" << ix << " jx::" << jx << " kx::" << kx << " lx::" << lx << std::endl;
      value += _int2e(ix,jx,kx,lx); if(m_type & 8) value -= _int2e(ix,kx,jx,lx);
      break;
    case 3:
//std::cout << "DEBUG[int2e_component]: 33" << std::endl;
      value += _int2e(ix,kx,jx,lx); if(m_type & 8) value -= _int2e(ix,jx,lx,kx);
      break;
    default:
      break;
  }
//std::cout << "DEBUG[int2e_component]: 34" << std::endl;
  // get scaling +/-
  return (m_type & 4) ? -value : value;
}
