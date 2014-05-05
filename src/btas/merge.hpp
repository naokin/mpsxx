

template<typename _T, unsigned long M, unsigned long N, class Q>
void merge(
const QSTArray<_T, M, Q>& A, const std::string& labelA,
      QSTArray<_T, N, Q>& B, const std::string& labelB)
{
}

void merge(
   const QSTArray<_T, M, Q>& A, const Index& indexA,
         QSTArray<_T, N, Q>& B, const Index& indexB)
{
   // this cannot be done at compile time to support string expression
   BTAS_RUNTIME_ASSERT(indexA.size() == M, "index size A must be the same as the rank of A");
   BTAS_RUNTIME_ASSERT(indexB.size() == N, "index size B must be the same as the rank of B");
}

merge(A, mergeInfo, B);

class MergeInfo {

public:


};

typedef unsigned long size_type;

// E.g.) A(i,j,k,l) -> B(ij,kl) : i,j,k,l = 0:3
// map info {ij}<-{i,j}
B.resize(info[0].size(), info[1].size());
for(size_type i = 0; i < info[0].size(); ++i) {
   if(info[0][i].extent() == 0) continue;

   for(size_type j = 0; j < info[1].size(); ++j) {
      if(info[1][j].extent() == 0) continue;

      B(i, j).resize(info[0][i].extent(), info[1][j].extent());
      B(i, j).fill(0.0);

      for(size_type ix = 0; ix < info[0][i].range().size(); ++ix) {
         Index ixm = info[0][i].index(ix);

         for(size_type jx = 0; jx < info[1][j].range().size(); ++jx) {
            Index jxm = info[1][j].index(jx);

            Slice<TArray<double, 2, Q>> b(B(i, j), info[0][i].range(ix), info[1][j].range(jx));
            b = A(ixm, jxm);
         }
      }
   }
}

typedef unsigned long size_type;

template<size_type _I, size_type _N, class = typename std::enable_if<(_I < _N)>::type>
struct __NDloop__Merge__0
{
   template<typename _T, size_type _M, size_type _N>
   __NDloop__Merge__0 (
      const Index<_N>& blocks,
      const Index<_N>& indexB,
      const Index<_N>& shapeB,
      const STArray<_T, _M>& A,
      const MergeInfo<_M, _N>& info,
            STArray<_T, _N>& B)
   {
      for (indexB[_I] = 0; indexB[_I] < blocks[_I]; ++indexB[_I]) {
         size_type n = info[_I][indexB[_I]].extent();
         if (n == 0) continue;
         shapeB[_I] = n;
         __NDloop__Merge__0<_I+1, _N> loop(blocks, indexB, shapeB, A, info, B);
      }
   }
};

template<size_type _N>
struct __NDloop__Merge__0<_N, _N>
{
   template<typename _T, size_type _M, size_type _N>
   __NDloop__Merge__0 (
      const Index<_N>& blocks,
      const Index<_N>& indexB,
      const Index<_N>& shapeB,
      const STArray<_T, _M>& A,
      const MergeInfo<_M, _N>& info,
            STArray<_T, _N>& B)
   {
      B(indexB).resize(shapeB, NumType<_T>::zero());
      Range<_N> rangeB = info.range(indexB);
      Index<_M> indexA;
      __NDloop__Merge__1<0, _N> loop(rangeB, indexA, A, B(indexB));
   }
};

struct __NDloop__Merge__1
{
   template<typename _T, size_type _M, size_type _N>
   __NDloop__Merge__0 (
      const Range<_N>& rangeB,
      const Index<_M>& indexA,
      const STArray<_T, _M>& A,
             TArray<_T, _N>& B)
   {
      for (size_type i = 0; i < rangeB[_I])
   }
};

#include <iostream>
#include <array>
#include <type_traits>

typedef unsigned long size_type;

template<size_type _I, size_type _N>
struct __NDloop {
   template<class Function, class = typename std::enable_if<(_I < _N)>::type>
   __NDloop (const std::array<size_type, _N>& range, std::array<size_type, _N>& index, Function fn) {
      for (index[_I] = 0; index[_I] < range[_I]; ++index[_I])
         __NDloop<_I+1, _N> loop(range, index, fn);
   }
};

template<size_type _N>
struct __NDloop<_N, _N> {
   template<class Function>
   __NDloop (const std::array<size_type, _N>& range, std::array<size_type, _N>& index, Function fn) { fn(index); }
};

template<size_type _N, class Function>
void NDloop (const std::array<size_type, _N>& range, Function fn) {
   std::array<size_type, _N> index;
   __NDloop<0, _N> loop(range, index, fn);
}

template<size_type _N>
void myFunc(const std::array<size_type, _N>& index) {
   std::cout << "[";
   for (size_type i = 0; i < _N; ++i) std::cout << " " << index[i];
   std::cout << " ]" << std::endl;
}

int main() {
   std::array<size_type, 4> range = { 2,2,2,2 };
   NDloop(range, myFunc<4>);
   return 0;
}
