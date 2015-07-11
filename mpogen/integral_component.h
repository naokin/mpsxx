#ifndef _MPSXX_CXX11_INTEGRAL_COMPONENT_H
#define _MPSXX_CXX11_INTEGRAL_COMPONENT_H 1

#include <legacy/DENSE/DArray.h>

#include "bit_operator_type.h"

namespace mpsxx     {

namespace fermionic {

/*! Determine 1-particle integral component for complementary operator construction in DMRG */
double int1e_component(const BIT_OPERATOR_TYPE& l_op, const BIT_OPERATOR_TYPE& s_op, const BIT_OPERATOR_TYPE& r_op, const btas::DArray<2>& _int1e);

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
double int2e_component(const BIT_OPERATOR_TYPE& l_op, const BIT_OPERATOR_TYPE& s_op, const BIT_OPERATOR_TYPE& r_op, const btas::DArray<4>& _int2e);

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_INTEGRAL_COMPONENT_H
