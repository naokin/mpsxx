#ifndef _UTILS_H
#define _UTILS_H 1

#include <cmath>

#include <btas/Dblas.h>

namespace util
{

template<int N>
void Normalize(btas::DArray<N>& x)
{
  double norm = btas::Ddot(x, x);
  btas::Dscal(1.0/sqrt(norm), x);
}

template<int N>
void Orthogonalize(const btas::DArray<N>& x, btas::DArray<N>& y)
{
  double ovlp = btas::Ddot(x, y);
  btas::Daxpy(-ovlp, x, y);
}

template<int N>
void Randomize(double noise, btas::DArray<N>& x, const btas::function<double(void)>& f_random_generator)
{
  btas::DArray<N> y(x.shape());
  y = f_random_generator;
  Normalize(y);
  btas::Daxpy(noise, y, x);
}

};

#endif // _UTILS_H
