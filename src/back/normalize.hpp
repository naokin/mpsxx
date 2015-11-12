#ifndef __MPSXX_DRIVER_NORMALIZE_HPP
#define __MPSXX_DRIVER_NORMALIZE_HPP

namespace mpsxx {

namespace driver {

template<class Q>
double Normalize(MPS<double, Q>& psi, bool forward = true)
{
   typedef typename MPS<double, Q>::size_type size_type;

   size_type n = psi.size();

   double nrm2;
   if(forward) {
      for(size_type i = 0; i < n-1; ++i) {
         btas:: STArray<double, 1>    s;
         btas::QSTArray<double, 3, Q> u;
         btas::QSTArray<double, 2, Q> v;
         btas::QSDgesvd(psi[i], u, s, v);
         psi[i].swap(u);

         btas::QSTArray<double, 3, Q> tmp;
         btas:: SDdidm(s, v);
         btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, v, psi[i+1], 1.0, tmp);
         psi[i+1].swap(tmp);
      }
      nrm2 = btas::QSDdotc(psi[n-1], psi[n-1]);
      btas::QSDscal(1.0/sqrt(nrm2), psi[n-1]);
   }
   else {
      for(size_type i = n-1; i > 0; --i) {
         btas:: STArray<double, 1>    s;
         btas::QSTArray<double, 2, Q> u;
         btas::QSTArray<double, 3, Q> v;
         btas::QSDgesvd(psi[i], u, s, v);
         psi[i].swap(v);

         btas::QSTArray<double, 3, Q> tmp;
         btas:: SDdimd(u, s);
         btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, psi[i-1], u, 1.0, tmp);
         psi[i-1].swap(tmp);
      }
      nrm2 = btas::QSDdotc(psi[0], psi[0]);
      btas::QSDscal(1.0/sqrt(nrm2), psi[0]);
   }
   return nrm2;
}

} // namespace driver

} // namespace mpsxx

#endif // __MPSXX_DRIVER_NORMALIZE_HPP
