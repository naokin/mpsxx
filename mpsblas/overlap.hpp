#ifndef __MPSXX_DMRG_OVERLAP_HPP
#define __MPSXX_DMRG_OVERLAP_HPP

#include "MPS_vector.hpp"
#include "MPO_vector.hpp"

namespace mpsxx {

template<typename T>
T overlap (const MPS_vector<T>& b, const MPS_vector<T>& a)
{
  using btas::shape;

  size_t N = b.size();
  assert(N == a.size());

  T ONE = static_cast<T>(1);

  btas::Tensor<T,4,CblasRowMajor> block;
  {
    btas::Tensor<T,4,CblasRowMajor> tmp4;
    btas::contract(ONE,b[0],shape(1),a[0],shape(1),ONE,block);
    btas::permute(tmp4,shape(0,2,1,3),block);
  }

  for(size_t i = 1; i < N; ++i) {
    btas::Tensor<T,5,CblasRowMajor> tmp5;
    btas::contract(ONE,block,shape(2),b[i],shape(0),ONE,tmp5);
    block.clear();
    btas::contract(ONE,tmp5,shape(2,3),a[i],shape(0,1),ONE,block);
  }

  assert(block.extent(0) == block.extent(2));
  assert(block.extent(1) == block.extent(3));
  T ovlp = static_cast<T>(0);
  for(size_t i = 0; i < block.extent(0); ++i)
    for(size_t j = 0; j < block.extent(1); ++j)
      ovlp += block(i,j,i,j);

  return ovlp;
}

} // namespace mpsxx

#endif // __MPSXX_DMRG_OVERLAP_HPP
