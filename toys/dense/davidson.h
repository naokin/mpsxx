#ifndef _DAVIDSON_H
#define _DAVIDSON_H 1

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <legacy/Dlapack.h>

#include "btas_template_specialize.h"
#include "utils.h"

namespace davidson
{

//
// Davidson's precondition
//
template<int N>
void precondition(double eval, const btas::DArray<N>& diag, btas::DArray<N>& errv)
{
  typename btas::DArray<N>::const_iterator id = diag.begin();
  for(typename btas::DArray<N>::iterator ir = errv.begin(); ir != errv.end(); ++ir, ++id) {
    double denm = eval - *id;
    if(fabs(denm) < 1.0e-12) denm = 1.0e-12;
    *ir /= denm;
  }
}

//
// Davidson eigen solver
//
template<int N>
double diagonalize(const boost::function<void(const btas::DArray<N>&, btas::DArray<N>&)>& f_contract,
                   const btas::DArray<N>& diag, btas::DArray<N>& wfnc)
{
  int max_ritz = 20;

  double eval = 0.0;

  // reserve working space
  std::vector< btas::DArray<N> > trial(max_ritz, btas::DArray<N>());
  std::vector< btas::DArray<N> > sigma(max_ritz, btas::DArray<N>());

  btas::Dcopy(wfnc, trial[0]);
  util::Normalize(trial[0]);
  f_contract(trial[0], sigma[0]);

  int niter = 0;
  int iconv = 0;
  while(iconv < 1 && niter < 20) {
    for(int m = 1; m <= max_ritz; ++m) {
      // compute small Hamiltonian matrix
      btas::DArray<2> heff(m, m);
      btas::DArray<2> ovlp(m, m);
      for(int i = 0; i < m; ++i) {
        heff(i, i) = btas::Ddot(trial[i], sigma[i]);
        ovlp(i, i) = btas::Ddot(trial[i], trial[i]);
        for(int j = 0; j < i; ++j) {
          double hij = btas::Ddot(trial[i], sigma[j]);
          heff(i, j) = hij;
          heff(j, i) = hij;
          double sij = btas::Ddot(trial[i], trial[j]);
          ovlp(i, j) = sij;
          ovlp(j, i) = sij;
        }
      }
      // solve eigenvalue problem to obtain Ritz value & vector
      btas::DArray<2> rvec;
      btas::DArray<1> rval;
      Dsyev(heff, rval, rvec);
      eval = rval(0);
      // rotate trial & sigma vectors by Ritz vector
      std::vector< btas::DArray<N> > trial_save(m, btas::DArray<N>());
      std::vector< btas::DArray<N> > sigma_save(m, btas::DArray<N>());
      for(int i = 0; i < m; ++i) {
        btas::Dcopy(trial[i], trial_save[i]);
        btas::Dcopy(sigma[i], sigma_save[i]);
        btas::Dscal(rvec(i, i), trial[i]);
        btas::Dscal(rvec(i, i), sigma[i]);
      }
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          if(i != j) {
            btas::Daxpy(rvec(i, j), trial_save[j], trial[i]);
            btas::Daxpy(rvec(i, j), sigma_save[j], sigma[i]);
          }
        }
      }
      // compute error vector
      btas::DArray<N> evec;
      btas::DArray<N> errv;
      btas::Dcopy( trial[0], evec);
      btas::Dcopy( sigma[0], errv);
      btas::Daxpy(-eval, evec, errv);
      double rnorm = btas::Ddot(errv, errv);
      if(rnorm < 1.0e-8) { ++iconv; break; }
      // solve correction equation
      if(m < max_ritz) {
        precondition(eval, diag, errv);
        for(int i = 0; i < m; ++i) {
          util::Normalize(errv);
          util::Orthogonalize(trial[i], errv);
        }
        util::Normalize(errv);
        btas::Dcopy(errv, trial[m]);
        sigma[m].free();
        f_contract(trial[m], sigma[m]);
      }
    }
    ++niter;
  }

  btas::Dcopy(trial[0], wfnc);

  return eval;
}

};

#endif // _DAVIDSON_H
