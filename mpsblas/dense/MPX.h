#ifndef __MPSXX_DMRG_DENSE_MPX_H
#define __MPSXX_DMRG_DENSE_MPX_H

#include <vector>

#include <legacy/DENSE/TArray.h>

namespace mpsxx {

/// MPO{0:Left,1:Bra,2:Ket,3:Right} -- Matrix Product Operators
template<typename T> using MPO = btas::TArray<T,4>;

template<typename T> using MPOs = std::vector<MPO<T>>;

/// MPS{0:Left,1:Ket,2:Right} -- Matrix Product States
template<typename T> using MPS = btas::TArray<T,3>;

template<typename T> using MPSs = std::vector<MPS<T>>;

/// Block{0:Bra,1:Op,2:Ket} -- renormalized BLOCK operator
template<typename T> using BLOCK = btas::TArray<T,3>;

} // namespace mpsxx

#endif // __MPSXX_DMRG_DENSE_MPX_H
