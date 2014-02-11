#ifndef __BTAS_REINDEX_HPP
#define __BTAS_REINDEX_HPP 1

#include <type_traits>

#include <btas_types.h>

namespace btas {

/// reindex (i.e. permute) for "any-rank" tensor
/// multiple loop is expanded at compile time
/// FIXME: how slower than explicit looping?
/// if considerably slower, should be specialized for small ranks (_N = 1 ~ 8?)
template<typename _T, size_type _N>
void Reindex (const _T* pX, _T* pY, const Index<_N>& strX, const Index<_N>& shapeY)
{
   __NDloop__reindex<0, _N> loop(pX, pY, 0, strX, shapeY);
}

/// NDloop class for Reindex
template<size_type _I, size_type _N, class = typename std::enable_if<(_I < _N-1)>::type>
struct __NDloop__reindex
{
   /// loop upon construction
   /// NOTE: pX and pY are passed as a reference of pointer to the next loop
   /// NOTE: on the other hand, addrX is passed as a value so that offset position (by addrX) is kept in this scope
   template<typename _T>
   __NDloop__reindex (const _T*& pX, _T*& pY, size_type addrX, const Index<_N>& strX, const Index<_N>& shapeY)
   {
      for (size_type i = 0; i < shapeY[_I]; ++i)
      {
         __NDloop__reindex<_I+1, _N> loop(pX, pY, addrX+i*strX[_I], strX, shapeY);
      }
   }
};

/// NDloop class for Reindex, specialized for the last index
template<size_type _N>
struct __NDloop_reindex<_N-1, _N>
{
   /// loop upon construction
   template<typename _T>
   __NDloop__reindex (const _T*& pX, _T*& pY, size_type addrX, const Index<_N>& strX, const Index<_N>& shapeY)
   {
      for (size_type i = 0; i < shapeY[_N-1]; ++i, ++pY)
      {
         *pY = pX[addrX+i*strX[_N-1]];
      }
   }
};

} // namespace btas

#endif // __BTAS_REINDEX_HPP
