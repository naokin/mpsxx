#ifndef _MPSXX_FERMIONIC_GET_PRODUCT_OPS_H
#define _MPSXX_FERMIONIC_GET_PRODUCT_OPS_H 1

#include <vector>

#include "bit_operator_type.h"

namespace mpsxx {

   namespace fermionic {

      std::vector<BIT_OPERATOR_TYPE> get_product_ops
         (const BIT_OPERATOR_TYPE& l_op, const BIT_OPERATOR_TYPE& s_op, const std::vector<size_t>& r_indxs, bool _swap_sweep_dir = false);

   };

};

#endif // _MPSXX_FERMIONIC_GET_PRODUCT_OPS_H
