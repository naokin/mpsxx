#include <iostream>

#include <cmath>
#include <algorithm>
#include <btas/Dblas.h>
#include <btas/Dlapack.h>
#include <btas/Dpermute.h>
#include "davidson.h"
#include "driver.h"
#include "driverlrt.h"
using namespace std;
using namespace btas;

void lrt::canonicalize         (bool forward,
                                const DArray<4>& sdm1,
                                const DArray<1>& val0,
                                const DArray<3>& mps0,
                                const DArray<3>& nul0,
                                      DArray<3>& mps1)
{
  mps1.resize(mps0.shape()); mps1 = 0.0;
  if(nul0.size() == 0) return;

  int M;
  if(forward) M = mps0.extent(2);
  else        M = mps0.extent(0);

  DArray<1> vali(M);
  for(int i = 0; i < M; ++i) {
    if(fabs(val0(i)) >= 1.0e-12) vali(i) = 1.0/val0(i);
    else                         vali(i) = 0.0;
  }

  if(forward) {
    DArray<3> scr1;
    Dgemm(NoTrans, NoTrans, 1.0, sdm1, mps0, 1.0, scr1);
    DArray<2> scr2;
    Dgemm(  Trans, NoTrans, 1.0, nul0, scr1, 1.0, scr2);
    Dgemm(NoTrans, NoTrans, 1.0, nul0, scr2, 1.0, mps1);
    Dleft_update(mps1, vali);
  }
  else {
    DArray<3> scr1;
    Dgemm(NoTrans,   Trans, 1.0, sdm1, mps0, 1.0, scr1);
    DArray<2> scr2;
    Dgemm(NoTrans, NoTrans, 1.0, nul0, scr1, 1.0, scr2);
    Dgemm(  Trans, NoTrans, 1.0, scr2, nul0, 1.0, mps1);
    Dright_update(vali, mps1);
  }
}

void lrt::compute_sigma_vector (const DArray<4>& mpo0,
                                const DArray<3>& wfn0, const DArray<3>& lstr0, const DArray<3>& rstr0,
                                const DArray<3>& wfn1, const DArray<3>& lstr1, const DArray<3>& rstr1,
                                      DArray<3>& sgv1)
{
  ::compute_sigma_vector(mpo0, lstr1, rstr0, wfn0, sgv1);
  ::compute_sigma_vector(mpo0, lstr0, rstr1, wfn0, sgv1);
  ::compute_sigma_vector(mpo0, lstr0, rstr0, wfn1, sgv1);
}

vector<double> lrt::diagonalize(const DArray<4>& mpo0,
                                const DArray<3>& diag,
                                const DArray<3>& lstr0,
                                const DArray<3>& rstr0,
                                const DArray<3>& wfn0,
                                const vector< DArray<3> >& lstr1,
                                const vector< DArray<3> >& rstr1,
                                      vector< DArray<3> >& wfn1,
                                int nroot, double tole, int max_iter)
{
  int max_ritz = max(20, nroot*2);
  // allocate working space
  vector<double> eigen_vals(nroot, 0.0);
  vector< DArray<3> > trial_vecs(max_ritz, DArray<3>(wfn0.shape()));
  vector< DArray<3> > sigma_vecs(max_ritz, DArray<3>(wfn0.shape()));
  vector< DArray<3> > lstr1_copy(max_ritz, DArray<3>(lstr0.shape()));
  vector< DArray<3> > rstr1_copy(max_ritz, DArray<3>(rstr0.shape()));
  // set left 1st-order block operators
  int nlstr1 = lstr1.size();
  for(int i = 0; i < nlstr1; ++i) {
    Dcopy(lstr1[i], lstr1_copy[i]);
  }
  for(int i = nlstr1; i < max_ritz; ++i) {
    lstr1_copy[i] = 0.0;
  }
  // set right 1st-order block operators
  int nrstr1 = rstr1.size();
  for(int i = 0; i < nrstr1; ++i) {
    Dcopy(rstr1[i], rstr1_copy[i]);
  }
  for(int i = nrstr1; i < max_ritz; ++i) {
    rstr1_copy[i] = 0.0;
  }
  // set 1st-order wavefunctions
  int nwfn1  = wfn1.size();
  for(int i = 0; i < nwfn1; ++i) {
    Dcopy(wfn1[i], trial_vecs[i]);
    Dorthogonalize(wfn0, trial_vecs[i]);
    Dnormalize(trial_vecs[i]);
    for(int j = 0; j < i; ++j) {
      Dorthogonalize(trial_vecs[j], trial_vecs[i]);
      Dnormalize(trial_vecs[i]);
    }
  }
  for(int i = nwfn1; i < nroot; ++i) {
    trial_vecs[i] = 0.0;
    compute_sigma_vector(mpo0, wfn0, lstr0, rstr0,
                         trial_vecs[i-1], lstr1_copy[i], rstr1_copy[i], trial_vecs[i]);
    Dorthogonalize(wfn0, trial_vecs[i]);
    Dnormalize(trial_vecs[i]);
    for(int j = 0; j < i; ++j) {
      Dorthogonalize(trial_vecs[j], trial_vecs[i]);
      Dnormalize(trial_vecs[i]);
    }
  }

  int iroot = 0;
  int iter  = 0;
  while(iter < max_iter && iroot < nroot) {
    for(int m = nroot; m <= max_ritz; ++m) {
      // compute sigma vectors up to nroot
      for(int i = 0; i < m; ++i) {
        sigma_vecs[i] = 0.0;
        compute_sigma_vector(mpo0, wfn0, lstr0, rstr0,
                             trial_vecs[i], lstr1_copy[i], rstr1_copy[i], sigma_vecs[i]);
      }
      // compute effective hamiltonian
      DArray<2> heff(m, m);
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          heff(i, j) = Ddot(trial_vecs[i], sigma_vecs[j]);
        }
      }
      // symmetrize
      for(int i = 0; i < m; ++i) {
        heff(i, i) *= 2.0;
        for(int j = 0; j < i; ++j) {
          heff(i, j) += heff(j, i);
          heff(j, i)  = heff(i, j);
        }
      }
      // diagonalize
      DArray<2> ritz_vec(m, m);
      DArray<1> ritz_val(m);
      Dsyev(heff, ritz_val, ritz_vec);
      for(int i = 0; i < nroot; ++i) {
        eigen_vals[i] = ritz_val(i);
      }
      // compute eigen vectors
      vector< DArray<3> > trial_vecs_save(m, DArray<3>());
      vector< DArray<3> > sigma_vecs_save(m, DArray<3>());
      for(int i = 0; i < m; ++i) {
        Dcopy(trial_vecs[i], trial_vecs_save[i]);
        Dcopy(sigma_vecs[i], sigma_vecs_save[i]);
        Dscal(ritz_vec(i, i), trial_vecs[i]);
        Dscal(ritz_vec(i, i), sigma_vecs[i]);
      }
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          if(i != j) {
            Daxpy(ritz_vec(i, j), trial_vecs_save[j], trial_vecs[i]);
            Daxpy(ritz_vec(i, j), sigma_vecs_save[j], sigma_vecs[i]);
          }
        }
      }
      // compute error vector
      DArray<3> eigvec;
      DArray<3> errvec;
      while(iroot < nroot) {
        Dcopy( trial_vecs[iroot], eigvec);
        Dcopy( sigma_vecs[iroot], errvec);
        Daxpy(-eigen_vals[iroot], eigvec, errvec);
        double rnorm = Ddot(errvec, errvec);
        if(rnorm < tole) break;
        ++iroot;
      }
      // solved
      if(iroot == nroot) break;
      // compute new trial vector
      if(m < max_ritz) {
        precondition(eigen_vals[iroot], diag, errvec);
        Dorthogonalize(wfn0, errvec);
        Dnormalize(errvec);
        for(int i = 0; i < m; ++i) {
          Dorthogonalize(trial_vecs[i], errvec);
          Dnormalize(errvec);
        }
        Dcopy(errvec, trial_vecs[m]);
      }
    }
    ++iter;
  }

  wfn1.resize(nroot, DArray<3>());
  for(int i = 0; i < nroot; ++i) {
    Dcopy(trial_vecs[i], wfn1[i]);
  }

  return eigen_vals;
}

//
// test for linear response
//

void lrt::compute_overlap_matrix(const vector< DArray<3> >& wfns,
                                 const vector< DArray<3> >& lmps,
                                 const vector< DArray<3> >& rmps)
{
  int L = wfns.size();
  Array< DArray<6>, 2 > S(L, L);
  for(int i = 0; i < L; ++i) {
    int li = wfns[i].extent(0);
    int ni = wfns[i].extent(1);
    int ri = wfns[i].extent(2);
    for(int j = 0; j < L; ++j) {
      int lj = wfns[j].extent(0);
      int nj = wfns[j].extent(1);
      int rj = wfns[j].extent(2);
      S(i, j).resize(li, ni, ri, lj, nj, rj);
      S(i, j) = 0.0;
      if(i == j) {
        for(int l = 0; l < li; ++l)
          for(int n = 0; n < ni; ++n)
            for(int r = 0; r < ri; ++r)
              S(i, j)(l, n, r, l, n, r) = 1.0;
      }
      if(i >  j) {
        DArray<6> I0;
        Dger(1.0, rmps[i], lmps[i-1], I0);
        DArray<6> I1;
        Dpermute(I0, shape(5, 1, 2, 3, 4, 0), I1);
        for(int k = i-1; k > j; --k) {
          DArray<6> I2;
          Dger(1.0, rmps[k], lmps[k-1], I2);
          DArray<6> I3;
          Dpermute(I2, shape(5, 1, 2, 3, 4, 0), I3);
          DArray<6> I4;
          Dgemm(NoTrans, NoTrans, 1.0, I1, I3, 1.0, I4);
          Dcopy(I4, I1);
        }
        Dcopy(I1, S(i, j));
      }
      if(i <  j) {
        DArray<6> I0;
        Dger(1.0, lmps[i], rmps[i+1], I0);
        DArray<6> I1;
        Dpermute(I0, shape(0, 1, 3, 2, 4, 5), I1);
        for(int k = i+1; k < j; ++k) {
          DArray<6> I2;
          Dger(1.0, lmps[k], rmps[k+1], I2);
          DArray<6> I3;
          Dpermute(I2, shape(0, 1, 3, 2, 4, 5), I3);
          DArray<6> I4;
          Dgemm(NoTrans, NoTrans, 1.0, I1, I3, 1.0, I4);
          Dcopy(I4, I1);
        }
        Dcopy(I1, S(i, j));
      }

      cout << "S[" << setw(2) << i << "," << setw(2) << j << "]: " << S(i, j).shape() << endl;
      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(2);
      for(int l0 = 0; l0 < li; ++l0) {
        for(int n0 = 0; n0 < ni; ++n0) {
          for(int r0 = 0; r0 < ri; ++r0) {
            for(int l1 = 0; l1 < lj; ++l1) {
              for(int n1 = 0; n1 < nj; ++n1) {
                for(int r1 = 0; r1 < rj; ++r1) {
                  cout << setw(6) << S(i, j)(l0, n0, r0, l1, n1, r1);
                }
              }
            }
            cout << endl;
          }
        }
      }
      cout << endl;
    }
  }
  cout << "contracted overlap matrix" << endl;
  DArray<2> S0(L, L);
  for(int i = 0; i < L; ++i) {
    for(int j = 0; j < L; ++j) {
      DArray<3> SC;
      Dgemm(NoTrans, NoTrans, 1.0, S(i, j), wfns[j], 1.0, SC);
      S0(i, j) = Ddot(wfns[i], SC);
      cout << setw(6) << S0(i, j);
    }
    cout << endl;
  }
  cout << endl;
}

