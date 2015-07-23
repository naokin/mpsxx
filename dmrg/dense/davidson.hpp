#ifndef __MPSXX_DMRG_DENSE_DAVIDSON_HPP
#define __MPSXX_DMRG_DENSE_DAVIDSON_HPP

#include <iostream>
#include <iomanip>

#include <vector>
#include <cmath>

#include <boost/function.hpp>

#include <legacy/DENSE/TArray.h>

namespace mpsxx {
namespace davidson {

/// Davidson's preconditioner
/// \param e approximate eigen value
/// \param d diagonal elements
/// \param r error vector
template<size_t N>
void precondition (double e, const btas::TArray<double,N>& d, btas::TArray<double,N>& r)
{
  auto di = d.begin();
  for(auto ri = r.begin(); ri != r.end(); ++ri, ++di) {
    double denom = e - *di;
    if(fabs(denom) < 1.0e-12)
      *ri  = 0.0;
    else
      *ri /= denom;
  }
}

/// Davidson eigen solver
/// \param fsigma function object to compute sigma vector
/// \param d diagonal elements
/// \param c on entry, guess vector. on exit, eigen vector
template<size_t N>
double diagonalize (
  const boost::function<void(
          const btas::TArray<double,N>&,
                btas::TArray<double,N>&)
        >& fsigma,
  const btas::TArray<double,N>& d,
        btas::TArray<double,N>& c)
{
  size_t MX_ITER = 4;
  size_t MX_RITZ = 20;

  double e = 0.0;

  // reserve working space
  std::vector<btas::TArray<double,N>> trial(MX_RITZ);
  std::vector<btas::TArray<double,N>> sigma(MX_RITZ);

  trial[0] = c;
  btas::Normalize(trial[0]);
  fsigma(trial[0],sigma[0]);

  size_t conv = 0;
  size_t iter = 0;

  size_t root = 1;

//std::cout.precision(8);
//std::cout.setf(std::ios::fixed,std::ios::floatfield);
//std::cout << "\t\t\t-------------------------------------------------------" << std::endl;
  while(conv < root && iter < MX_ITER) {
    for(size_t m = 1; m <= MX_RITZ; ++m) {
      // compute small Hamiltonian matrix
      btas::TArray<double,2> H(m, m);
      btas::TArray<double,2> S(m, m);
      for(size_t i = 0; i < m; ++i) {
        H(i, i) = btas::Dot(trial[i], sigma[i]);
        S(i, i) = btas::Dot(trial[i], trial[i]);
        for(size_t j = 0; j < i; ++j) {
          double Hij = btas::Dot(trial[i],sigma[j]);
          H(i, j) = Hij;
          H(j, i) = Hij;
          double Sij = btas::Dot(trial[i],trial[j]);
          S(i, j) = Sij;
          S(j, i) = Sij;
        }
      }
    //std::cout << "DEBUG :: H = " << H << std::endl;
    //std::cout << "DEBUG :: S = " << S << std::endl;
      // solve eigenvalue problem to obtain Ritz value & vector
      btas::TArray<double,2> rvec;
      btas::TArray<double,1> rval;
      btas::Syev('V','U',H,rval,rvec);
      e = rval(0);

      // rotate trial & sigma vectors by Ritz vector
      std::vector<btas::TArray<double,N>> trisave(m);
      std::vector<btas::TArray<double,N>> sgvsave(m);
      for(size_t i = 0; i < m; ++i) {
        trisave[i] = trial[i];
        sgvsave[i] = sigma[i];
        btas::Scal(rvec(i,i),trial[i]);
        btas::Scal(rvec(i,i),sigma[i]);
      }
      for(size_t i = 0; i < m; ++i)
        for(size_t j = 0; j < m; ++j)
          if(i != j) {
//          btas::Axpy(rvec(i,j),trisave[j],trial[i]);
//          btas::Axpy(rvec(i,j),sgvsave[j],sigma[i]);
            btas::Axpy(rvec(j,i),trisave[j],trial[i]);
            btas::Axpy(rvec(j,i),sgvsave[j],sigma[i]);
          }

      // compute error vector
      btas::TArray<double,N> v(trial[0]);
      btas::TArray<double,N> r(sigma[0]);
      btas::Axpy(-e,v,r);

      double rnorm = btas::Dot(r,r);
    //std::cout << "\t\t\tmicro iteractions [" << std::setw(2) << m << "] :: " << std::setw(12) << e << " ( " << std::setw(12) << rnorm << " ) " << std::endl;
      if(rnorm < 1.0e-8) {
        ++conv; break;
      }
      else {
        // solve correction equation
        if(m < MX_RITZ) {
          precondition(e,d,r);
          for(size_t i = 0; i < m; ++i) {
            Normalize(r);
            Orthogonalize(trial[i],r);
          }
          Normalize(r);
          trial[m] = r;
          sigma[m].clear();
          fsigma(trial[m],sigma[m]);
        }
      }
    }
    ++iter;
  }
//std::cout << "\t\t\t-------------------------------------------------------" << std::endl;
  c = trial[0];
  return e;
}

} // namespace davidson
} // namespace mpsxx

#endif // __MPSXX_DMRG_DENSE_DAVIDSON_HPP
