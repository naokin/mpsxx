#ifndef __MPSXX_TYPEDEFS_H
#define __MPSXX_TYPEDEFS_H

#include <vector>
#include <MP/MPS.hpp>
#include <MP/MPO.hpp>

namespace mpsxx {

template<class Q, class Matrix>
using MPS_vector = std::vector<MPS<Q,Matrix>>;

template<class Q, class Matrix>
using MPO_vector = std::vector<MPO<Q,Matrix>>;

} // namespace mpsxx

#endif // __MPSXX_TYPEDEFS_H
