#ifndef PROTOTYPE_DAVIDSON_H
#define PROTOTYPE_DAVIDSON_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <btas/Dblas.h>
#include <btas/Dlapack.h>

extern double random_gen(void);

// davidson's precondition
template<int N>
void precondition(double e, const btas::DArray<N>& diag, btas::DArray<N>& errv)
{
  typename btas::DArray<N>::const_iterator idiag = diag.begin();
  typename btas::DArray<N>::iterator       ierrv = errv.begin();
  for(; ierrv != errv.end(); ++idiag, ++ierrv) {
    double denm = e - *idiag;
    if(fabs(denm) >= 1.0e-8) *ierrv /= denm;
  }
}

// olsen's precondition
template<int N>
void precondition(double e, const btas::DArray<N>& diag, btas::DArray<N>& errv, const btas::DArray<N>& eigv)
{
  btas::DArray<N> eig0;
  btas::Dcopy(eigv, eig0);
  precondition(e, diag, eig0);
  double numr = btas::Ddot(eig0, errv);
  double denm = btas::Ddot(eigv, eig0);
  btas::Daxpy(-numr/denm, eigv, errv);
  precondition(e, diag, errv);
}

// davidson eigen solver
template<int N>
std::vector<double> davidson(const boost::function<void(const btas::DArray<N>&, btas::DArray<N>&)>& sgv_functor,
                             const btas::DArray<N>& diag, std::vector< btas::DArray<N> >& eigen_vecs,
                             int nroot = 1, double tole = 1.0e-8)
{
  int max_iter = 100;
  int max_ritz = std::max(20, nroot*2);
  if(max_ritz > diag.size()) max_ritz = diag.size();

  std::vector<double> eigen_vals(nroot, 0.0);

  const btas::TinyVector<int, N>& psi_shape = diag.shape();
  std::vector< btas::DArray<N> > trial_vecs(max_ritz, btas::DArray<N>(psi_shape));
  std::vector< btas::DArray<N> > sigma_vecs(max_ritz, btas::DArray<N>(psi_shape));

  for(int i = 0; i < eigen_vecs.size(); ++i) {
    btas::Dcopy(eigen_vecs[i], trial_vecs[i]);
    for(int j = 0; j < i; ++j) {
      btas::Dorthogonalize(trial_vecs[j], trial_vecs[i]);
    }
    btas::Dnormalize(trial_vecs[i]);
  }
  for(int i = eigen_vecs.size(); i < nroot; ++i) {
    trial_vecs[i] = random_gen;
    for(int j = 0; j < i; ++j) {
      btas::Dorthogonalize(trial_vecs[j], trial_vecs[i]);
    }
    btas::Dnormalize(trial_vecs[i]);
  }

  for(int m = 0; m < nroot; ++m) {
    sigma_vecs[m] = 0.0;
    sgv_functor(trial_vecs[m], sigma_vecs[m]);
  }

  int iroot = 0;
  int iter  = 0;
  while(iter < max_iter && iroot < nroot) {
    for(int m = nroot; m <= max_ritz; ++m) {
      btas::DArray<2> heff(m, m);
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          heff(j, i) = btas::Ddot(trial_vecs[j], sigma_vecs[i]);
        }
      }

      btas::DArray<2> ritz_vec(m, m);
      btas::DArray<1> ritz_val(m);
      btas::Dsyev(heff, ritz_val, ritz_vec);
      for(int i = 0; i < nroot; ++i) {
        eigen_vals[i] = ritz_val(i);
      }

      std::vector< btas::DArray<N> > trial_vecs_save(m, btas::DArray<N>());
      std::vector< btas::DArray<N> > sigma_vecs_save(m, btas::DArray<N>());
      for(int i = 0; i < m; ++i) {
        btas::Dcopy(trial_vecs[i], trial_vecs_save[i]);
        btas::Dcopy(sigma_vecs[i], sigma_vecs_save[i]);
        btas::Dscal(ritz_vec(i, i), trial_vecs[i]);
        btas::Dscal(ritz_vec(i, i), sigma_vecs[i]);
      }
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          if(i != j) {
            btas::Daxpy(ritz_vec(i, j), trial_vecs_save[j], trial_vecs[i]);
            btas::Daxpy(ritz_vec(i, j), sigma_vecs_save[j], sigma_vecs[i]);
          }
        }
      }

      btas::DArray<N> eigvec;
      btas::DArray<N> errvec;
      while(iroot < nroot) {
        btas::Dcopy( trial_vecs[iroot], eigvec);
        btas::Dcopy( sigma_vecs[iroot], errvec);
        btas::Daxpy(-eigen_vals[iroot], eigvec, errvec);
        double rnorm = btas::Ddot(errvec, errvec);
        if(rnorm >= tole) break;
        ++iroot;
      }
      if(iroot == nroot) break;
      if(m < max_ritz) {
        precondition(eigen_vals[iroot], diag, errvec);
        for(int i = 0; i < m; ++i) {
          btas::Dnormalize(errvec);
          btas::Dorthogonalize(trial_vecs[i], errvec);
        }
        btas::Dnormalize(errvec);
        btas::Dcopy(errvec, trial_vecs[m]);

        sigma_vecs[m] = 0.0;
        sgv_functor(trial_vecs[m], sigma_vecs[m]);
      }
    }
    ++iter;
  }

//std::cout << "\t\t\tDEBUG::davidson converged after " << std::setw(3) << iter << " iterations" << std::endl;

  eigen_vecs.resize(nroot, btas::DArray<N>(psi_shape));
  for(int i = 0; i < nroot; ++i) {
    btas::Dcopy(trial_vecs[i], eigen_vecs[i]);
  }

  return eigen_vals;
}

// davidson eigen solver for lowest eigenvalue
template<int N>
double davidson(const boost::function<void(const btas::DArray<N>&, btas::DArray<N>&)>& sgv_functor,
                const btas::DArray<N>& diag, btas::DArray<N>& eigen_vec, double tole = 1.0e-8)
{
  std::vector< btas::DArray<N> > eigen_vecs(1, btas::DArray<N>());
  btas::Dcopy(eigen_vec, eigen_vecs[0]);
  std::vector<double> eigen_vals;
  eigen_vals = davidson(sgv_functor, diag, eigen_vecs, 1, tole);
//eigen_vals = davidson(sgv_functor, diag, eigen_vecs, 4, tole);
//std::cout << "\t\t\tDEBUG::davidson printing eigen values" << std::endl;
//std::cout.setf(std::ios::fixed, std::ios::floatfield);
//std::cout.precision(16);
//for(int i = 0; i < eigen_vals.size(); ++i) {
//  std::cout << "\t\t\t\teigen value[" << std::setw(2) << i << "]: " << std::setw(24) << eigen_vals[i] << std::endl;
//}
  btas::Dcopy(eigen_vecs[0], eigen_vec);
  return eigen_vals[0];
}

#endif // PROTOTYPE_DAVIDSON_H
