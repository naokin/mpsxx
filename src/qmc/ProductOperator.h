#ifndef __MPSXX_QMC_PRODUCT_OPERATOR_H
#define __MPSXX_QMC_PRODUCT_OPERATOR_H 1

namespace mpsxx {

template<class Q, unsigned long N>
struct ProductOperator {
   double coeff_;
   btas::IVector<N> index_;
   btas::TVector<btas::QSTArray<double, 2, Q>, N> ops_;
};

template<class Q>
struct ProductOperator<Q, 0> {
   double coeff_;
};

template<unsigned long I, unsigned long N, class = typename std::enable_if<(I < N)>::type>
struct __loop_over_operators
{
   template<class Q>
   __loop_over_operators (MPS<Q>& x, const ProductOperator<Q, N>& op)
   {
      btas::
      __loop_over_operator<I+1, N> loop(x, op);
   }
};

template<unsigned long N>
struct __loop_over_operators<N, N>
{
   template<class Q>
   __loop_over_operators (MPS<Q>& x, const ProductOperator<Q, N>& op) { }
};

} // namespace mpsxx

#endif // __MPSXX_QMC_PRODUCT_OPERATOR_H
