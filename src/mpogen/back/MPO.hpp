#ifndef __MPSXX_MPO_MPO_HPP
#define __MPSXX_MPO_MPO_HPP

#include <vector>

namespace mpsxx {

#ifdef _ENABLE_QUANTA
template<typename Tp, class Q>
using MPO = std::vector<btas::QSTArray<Tp, 4, Q>>;
#else
template<typename Tp>
using MPO = std::vector<btas::TArray<Tp, 4>>;
#endif // _ENABLE_QUANTA

} // namespace mpsxx

#endif // __MPSXX_MPO_MPO_HPP
