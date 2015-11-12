#ifndef __MPSXX_DMRG_MPS_VECTOR_HPP
#define __MPSXX_DMRG_MPS_VECTOR_HPP

#include <vector>
#include <btas/tensor.hpp>

namespace mpsxx {

/// rank-3 tensor represents A{al,ni,ar}
template<typename T>
using MPS = btas::Tensor<T,3,CblasRowMajor>;

template<typename T>
using MPS_vector = std::vector<MPS<T>>;

} // namespace mpsxx

#endif // __MPSXX_DMRG_MPS_VECTOR_HPP
