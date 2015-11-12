#ifndef __MPSXX_DMRG_EXPECTATION_HPP
#define __MPSXX_DMRG_EXPECTATION_HPP

#include "MPS_vector.hpp"
#include "MPO_vector.hpp"

namespace mpsxx {

template<typename T>
T expectation (const MPS_vector<T>& b, const MPO_vector<T>& H, const MPS_vector<T>& a)
{
  using btas::shape;

  size_t N = H.size();
  assert(N == b.size());
  assert(N == a.size());

  T ONE = static_cast<T>(1);

  btas::Tensor<T,6,CblasRowMajor> block;
  {
    btas::Tensor<T,5,CblasRowMajor> tmp5;
    btas::contract(ONE,b[0],shape(1),H[0],shape(1),ONE,tmp5);
    btas::Tensor<T,6,CblasRowMajor> tmp6;
    btas::contract(ONE,tmp5,shape(3),a[0],shape(1),ONE,tmp6);
    btas::permute(tmp6,shape(0,2,4,1,3,5),block);
  }

  for(size_t i = 1; i < N; ++i) {
    btas::Tensor<T,7,CblasRowMajor> tmp7_1;
    btas::contract(ONE,block,shape(3),b[i],shape(0),ONE,tmp7_1);
    btas::Tensor<T,7,CblasRowMajor> tmp7_2;
    btas::contract(ONE,tmp7_1,shape(3,5),H[i],shape(0,1),ONE,tmp7_2);
    block.clear();
    btas::contract(ONE,tmp7_2,shape(3,5),a[i],shape(0,1),ONE,block);
  }

  assert(block.extent(0) == block.extent(3));
  assert(block.extent(1) == block.extent(4));
  assert(block.extent(2) == block.extent(5));
  T expt = static_cast<T>(0);
  for(size_t i = 0; i < block.extent(0); ++i)
    for(size_t j = 0; j < block.extent(1); ++j)
      for(size_t k = 0; k < block.extent(2); ++k)
        expt += block(i,j,k,i,j,k);

  return expt;
}

} // namespace mpsxx

#endif // __MPSXX_DMRG_EXPECTATION_HPP
