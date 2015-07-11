#ifndef __MPSXX_MPS_H
#define __MPSXX_MPS_H 1

#include <vector>

#include <legacy/QSTArray.h>

namespace mpsxx {

template<class Q>
using MPS = std::vector<btas::QSTArray<double, 3, Q>>;

template<class Q>
using MPO = std::vector<btas::QSTArray<double, 4, Q>>;

template<class Q>
double dot (const MPS<Q>& x, const MPS<Q>& y)
{
   typedef typename MPS<Q>::size_type size_type;

   const size_type L = x.size();
   MPSXX_RUNTIME_ASSERT(L == y.size(), "number of sites must be the same");

   btas::Qshape<Q> q0(1, Q::zero());
   btas::Dshape    d0(q0.size(), 1);
   btas::QSTArray<double, 2, Q> left_(Q::zero(), btas::make_array(q0,q0), btas::make_array(d0,d0), 1.0);

   for(size_type i = 0; i < L; ++i) {
      btas::QSTArray<double, 2, Q> f_;
      btas::QSTgemm(btas::NoTrans, btas::NoTrans, 1.0, left_, x[i], 1.0, f_);
      left_.clear();
      btas::QSTgemm(btas::ConjTrans, btas::NoTrans, 1.0, y[i], f_, 1.0, left_);
   }

   return left_(0,0)(0,0);
}

template<class Q>
double dot (const MPS<Q>& x, const MPO<Q>& h, const MPS<Q>& y)
{
   typedef typename MPS<Q>::size_type size_type;

   const size_type L = x.size();
   MPSXX_RUNTIME_ASSERT(L == h.size(), "number of sites must be the same");
   MPSXX_RUNTIME_ASSERT(L == y.size(), "number of sites must be the same");

   btas::Qshape<Q> q0(1, Q::zero());
   btas::Dshape    d0(q0.size(), 1);
   btas::QSTArray<double, 3, Q> left_(Q::zero(), btas::make_array(q0,q0,q0), btas::make_array(d0,d0,d0), 1.0);

   for(size_type i = 0; i < L; ++i) {
      btas::QSTArray<double, 4, Q> g_;
      btas::QSTgemm(btas::ConjTrans, btas::NoTrans, 1.0, x[i], left_, 1.0, g_);

      btas::QSTArray<double, 4, Q> f_;
      btas::QSTcontract(1.0, g_, btas::shape(2,0), h, btas::shape(0,1), 1.0, f_);

      left_.clear();
      btas::QSTcontract(1.0, f_, btas::shape(1,2), y, btas::shape(0,1), 1.0, left_);
   }

   return left_(0,0,0)(0,0,0);
}

template<class Q>
void scal (const double& alpha, MPS<Q>& x)
{
   if (x.size() == 0) return;
   btas::QSTscal(alpha, x[0]);
}

} // namespace mpsxx

#endif // __MPSXX_MPS_H
