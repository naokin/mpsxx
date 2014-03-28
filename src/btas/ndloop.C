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
