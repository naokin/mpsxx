#ifndef __MPSXX_MPX_TYPES_H
#define __MPSXX_MPX_TYPES_H

#include <vector>

#include <legacy/QSPARSE/QSTArray.h>

namespace mpsxx {

/// MPO{0:Left,1:Bra,2:Ket,3:Right} -- Matrix Product Operators
template<typename T, class Q> using MPO = btas::QSTArray<T,4,Q>;

template<typename T, class Q> using MPOs = std::vector<MPO<T,Q>>;

/// MPS{0:Left,1:Ket,2:Right} -- Matrix Product States
template<typename T, class Q> using MPS = btas::QSTArray<T,3,Q>;

template<typename T, class Q> using MPSs = std::vector<MPS<T,Q>>;

/// Block{0:Bra,1:Op,2:Ket} -- renormalized BLOCK operator
template<typename T, class Q> using BLOCK = btas::QSTArray<T,3,Q>;

} // namespace mpsxx

#endif // __MPSXX_MPO_H
