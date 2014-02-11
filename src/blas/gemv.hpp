#ifndef __MPSXX_BLAS_GEMV_HPP
#define __MPSXX_BLAS_GEMV_HPP 1

namespace mpsxx {

namespace blas {

template<class Q>
void gemv(const double& alpha, const MPO<Q>& A, const MPS<Q>& X, const double& beta, MPS<Q>& Y)
{
   typedef typename MPS<Q>::size_type size_type;

   size_type L = X.size(); MPSXX_RUNTIME_ASSERT(A.size() == L, "number of sites (A vs. X) must be the same");

   if(Y.size() == 0) {
      Y.resize(L);
   }
   else {
      // test size of Y is equal to those of A and X
      MPSXX_RUNTIME_ASSERT(Y.size() == L, "number of sites (Y) must be the same");
      // test whether components of Y are either "all empty" or "all allocated"
      bool is_all_empty = true;
      bool is_all_alloc = true;
      for(size_type i = 0; i < L; ++i) {
         bool is_empty = Y[i].empty();
         is_all_empty &=  is_empty;
         is_all_alloc &= !is_empty;
      }
      MPSXX_RUNTIME_ASSERT(is_all_empty || is_all_alloc, "output MPS has missing component");
   }

   // Ti is used as a "bare contracted" tensor
   btas::QSTArray<double, 5, Q> Ti;
   // alpha is multiplied at 0-th component of AX
   btas::QSTIndexedContract(alpha, A[0], "bl,px,qx,br", X[0], "al,qx,ar", 1.0, Ti, "bl,al,px,br,ar");
   // beta is multiplied at 0-th component of Y
   if(!Y[0].empty()) {
      btas::QSTscal(beta, Y[0]);
   }

   // Tm is used as a "half-merged" tensor
   btas::QSTArray<double, 4, Q> Tm;
   // merge QN of left-edge (OBC is supposed)
   btas::QSTIndexedMerge(Ti, "bl,al,px,br,ar", Tm, "{bl,al},px,br,ar");

   for(size_type i = 1; i < L; ++i) {
      Ti.clear();
      // MPO x MPS
      btas::QSTIndexedContract(1.0, A[i], "bl,px,qx,br", X[i], "al,qx,ar", 1.0, Ti, "bl,al,px,br,ar");

      // take copy of Tm
      btas::QSTArray<double, 4, Q> AX(Tm);
      btas::QSTArray<double, 3, Q> Yi;
      // merge QN simultaneously with two tensors to form two merged tensors
      btas::QSTIndexedMerge(AX, "al,px,bm,am", Ti, "bm,am,qx,br,ar", Yi, "al,px,{bm,am}", Tm, "{bm,am},qx,br,ar");

      //
      if(Y[i].empty()) {
         Y[i-1] = Yi;
      }
      else {
         btas::QSTIndexedDsum(Y[i-1], "al,px,ar", Yi, "bl,px,br", Yx);
      }
   }
}

} // namespace blas

} // namespace mpsxx

#endif // __MPSXX_BLAS_GEMV_HPP
