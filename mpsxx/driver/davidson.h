#ifndef _PROTOTYPE_DAVIDSON_H
#define _PROTOTYPE_DAVIDSON_H 1

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <legacy/QSPARSE/QSTArray.h>

namespace davidson {

template<size_t N, class Q>
using Functor = boost::function<void(const btas::QSTArray<double, N, Q>&, btas::QSTArray<double, N, Q>&)>;

//
// Davidson's precondition
//

template<size_t N, class Q>
void precondition
(const double& eval, const btas::QSTArray<double, N, Q>& diag, btas::QSTArray<double, N, Q>& errv)
{
  for(auto ir = errv.begin(); ir != errv.end(); ++ir) {
    auto id = diag.find(ir->first);
    if(id != diag.end()) {
      auto irx = ir->second->begin();
      auto idx = id->second->begin();
      for(; irx != ir->second->end(); ++irx, ++idx) {
        double denm = eval - *idx;
        if(fabs(denm) < 1.0e-12) denm = 1.0e-12;
        *irx /= denm;
      }
    }
    else {
      btas::Dscal(1.0/eval, (*ir->second));
    }
  }
}

//
// Davidson eigen solver
//

template<size_t N, class Q>
double diagonalize
(const Functor<N, Q>& f_contract, const btas::QSTArray<double, N, Q>& diag, btas::QSTArray<double, N, Q>& wfnc)
{
  int max_ritz = 10;

  double eval = 0.0;

  // reserve working space
  std::vector<btas::QSTArray<double, N, Q>> trial(max_ritz);
  std::vector<btas::QSTArray<double, N, Q>> sigma(max_ritz);

  btas::Copy(wfnc, trial[0]);

  int niter = 0;
  int iconv = 0;
  while(iconv < 1 && niter < 4) {
    std::cout << "\t\t\tmacro iteration [ " << std::setw(2) << niter << " ] " << std::endl;
    std::cout << "\t\t\t--------------------------------------------------" << std::endl;
    // to keep numerical stability
    btas::Normalize(trial[0]);
    sigma[0].clear();
    f_contract(trial[0], sigma[0]);
    for(int m = 1; m <= max_ritz; ++m) {
      // compute small Hamiltonian matrix
      btas::TArray<double, 2> heff(m, m);
      btas::TArray<double, 2> ovlp(m, m);
      for(int i = 0; i < m; ++i) {
        heff(i, i) = btas::Dotc(trial[i], sigma[i]);
        ovlp(i, i) = btas::Dotc(trial[i], trial[i]);
        for(int j = 0; j < i; ++j) {
          double hij = btas::Dotc(trial[i], sigma[j]);
          heff(i, j) = hij;
          heff(j, i) = hij;
          double sij = btas::Dotc(trial[i], trial[j]);
          ovlp(i, j) = sij;
          ovlp(j, i) = sij;
        }
      }
      // solve eigenvalue problem to obtain Ritz value & vector
      btas::TArray<double, 2> rvec;
      btas::TArray<double, 1> rval;
      Syev('V', 'U', heff, rval, rvec);
      eval = rval(0);
      std::cout << "\t\t\tmicro iteration [ " << std::setw(2) << m << " ] :: " << std::setprecision(16) << std::setw(24) << std::fixed << eval << std::endl;
      // rotate trial & sigma vectors by Ritz vector
      std::vector<btas::QSTArray<double, N, Q>> trial_save(m);
      std::vector<btas::QSTArray<double, N, Q>> sigma_save(m);
      for(int i = 0; i < m; ++i) {
        btas::Copy(trial[i], trial_save[i]);
        btas::Copy(sigma[i], sigma_save[i]);
        btas::Scal(rvec(i, i), trial[i]);
        btas::Scal(rvec(i, i), sigma[i]);
      }
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          if(i != j) {
            btas::Axpy(rvec(i, j), trial_save[i], trial[j]);
            btas::Axpy(rvec(i, j), sigma_save[i], sigma[j]);
          }
        }
      }
      // compute error vector
      btas::QSTArray<double, N, Q> evec;
      btas::QSTArray<double, N, Q> errv;
      btas::Copy( trial[0], evec);
      btas::Copy( sigma[0], errv);
      btas::Axpy(-eval, evec, errv);
      double rnorm = btas::Dotc(errv, errv);
      if(rnorm < 1.0e-8) { ++iconv; break; }
      // solve correction equation
      if(m < max_ritz) {
        precondition(eval, diag, errv);
        for(int i = 0; i < m; ++i) {
          btas::Normalize(errv);
          btas::Orthogonalize(trial[i], errv);
        }
        btas::Normalize(errv);
        btas::Copy(errv, trial[m]);
        sigma[m].clear();
        f_contract(trial[m], sigma[m]);
      }
    }
    ++niter;
    std::cout << "\t\t\t--------------------------------------------------" << std::endl;
  }
  btas::Copy(trial[0], wfnc);

  return eval;
}

};

#endif // _PROTOTYPE_DAVIDSON_H
