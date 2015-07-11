#ifndef _MPSXX_CXX11_GENERATE_SITE_OPERATOR_H
#define _MPSXX_CXX11_GENERATE_SITE_OPERATOR_H 1

#include <legacy/DENSE/DArray.h>
#include <legacy/QSPARSE/QSDArray.h>

#include "bit_operator_type.h"

namespace mpsxx     {

   namespace fermionic {

      /*!
       *  Generate matrix product operator on single site
       */
      void generate_site_operator
         (btas::QSDArray<4, Quantum>& op,
          const BIT_OPERATOR_TYPE& l_op, const size_t& l_index,
          const BIT_OPERATOR_TYPE& s_op,
          const BIT_OPERATOR_TYPE& r_op, const size_t& r_index,
          const btas::DArray<2>& oneint, const btas::DArray<4>& twoint);

   }; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_GENERATE_SITE_OP_H
