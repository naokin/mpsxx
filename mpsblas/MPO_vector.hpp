#ifndef __MPSXX_DMRG_MPO_VECTOR_HPP
#define __MPSXX_DMRG_MPO_VECTOR_HPP

#include <vector>
#include <btas/tensor.hpp>

namespace mpsxx {

/// rank-4 tensor represents W{bl,ni',ni,br}
template<typename T>
using MPO = btas::Tensor<T,4,CblasRowMajor>;

template<typename T>
using MPO_vector = std::vector<MPO<T>>;

} // namespace mpsxx

#endif // __MPSXX_DMRG_MPO_VECTOR_HPP
