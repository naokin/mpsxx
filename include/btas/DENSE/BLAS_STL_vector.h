#ifndef __BTAS_BLAS_STL_VECTOR_H
#define __BTAS_BLAS_STL_VECTOR_H 1

#include <vector>

#include <btas/common/btas.h> // BTAS_THROW

#include <blas/wrappers.h>

namespace btas
{

template<typename T>
void Copy (const std::vector<T>& x, std::vector<T>& y)
{
   y.resize(x.size());

   copy(x.size(), x.data(), 1, y.data(), 1);
}

template<typename T>
void Scal (const T& alpha, std::vector<T>& x)
{
   scal(x.size(), alpha, x.data(), 1);
}

template<typename T>
void Axpy (const T& alpha, const std::vector<T>& x, std::vector<T>& y)
{
   if(y.size() > 0)
   {
      BTAS_THROW(x.size() == y.size(), "Axpy(std::vector): x and y must have the same size.");
   }
   else
   {
      y.resize(x.size(), static_cast<T>(0));
   }

   axpy(x.size(), alpha, x.data(), 1, y.data(), 1);
}

} // namespace btas

#endif // __BTAS_BLAS_STL_VECTOR_H
