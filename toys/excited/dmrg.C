#include <iostream>
#include <iomanip>
#include <algorithm>
#include <driver/davidson.h>
#include <btas/Dpermute.h>
#include <btas/Dcontract.h>
#include "driver.h"
#include "dmrg.h"
using namespace std;
using namespace btas;

extern double random_gen(void);

void mpo_init(const DmrgInput& input, MpStorages& sites)
{
  int L = input.N_sites;
  int d = input.N_phys_index;
  double J  = input.J;
  double Jz = input.Jz;
  double hz = input.hz;
  double sz = static_cast<double>(d - 1) / 2;
  // resize & set to 0
  sites[0].mpopr.resize(1, d, d, 5);
  sites[0].mpopr = 0;
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpopr.resize(5, d, d, 5);
    sites[i].mpopr = 0;
  }
  sites[L-1].mpopr.resize(5, d, d, 1);
  sites[L-1].mpopr = 0;
  // construct mpos for spin-hamiltonian
  double mz = sz;
  for(int m = 0; m < d; ++m) {
    sites[0].mpopr(0, m,   m,   0) =  hz * mz; // hz  Sz
    sites[0].mpopr(0, m,   m,   3) =  Jz * mz; // Jz  Sz
    sites[0].mpopr(0, m,   m,   4) =  1.0;     // I
    for(int i = 1; i < L-1; ++i) {
      sites[i].mpopr(0, m,   m,   0) =  1.0;     //     I
      sites[i].mpopr(3, m,   m,   0) =  mz;      //     Sz
      sites[i].mpopr(4, m,   m,   0) =  hz * mz; // hz  Sz
      sites[i].mpopr(4, m,   m,   3) =  Jz * mz; // Jz  Sz
      sites[i].mpopr(4, m,   m,   4) =  1.0;     // I
    }
    sites[L-1].mpopr(0, m,   m,   0) =  1.0;     //     I
    sites[L-1].mpopr(3, m,   m,   0) =  mz;      //     Sz
    sites[L-1].mpopr(4, m,   m,   0) =  hz * mz; // hz  Sz
    mz -= 1.0;
  }
  double mz_plus  = sz - 1.0;
  double mz_minus = sz;
  for(int m = 0; m < d - 1; ++m) {
    double c_plus  = sqrt(sz*(sz+1.0)-mz_plus *(mz_plus +1.0));
    double c_minus = sqrt(sz*(sz+1.0)-mz_minus*(mz_minus-1.0));
    sites[0].mpopr(0, m,   m+1, 1) =  c_plus  * J / 2;   // J/2 S+
    sites[0].mpopr(0, m+1, m,   2) =  c_minus * J / 2;   // J/2 S-
    for(int i = 1; i < L-1; ++i) {
      sites[i].mpopr(1, m+1, m,   0) =  c_minus * 1.0;     //     S-
      sites[i].mpopr(2, m,   m+1, 0) =  c_plus  * 1.0;     //     S+
      sites[i].mpopr(4, m,   m+1, 1) =  c_plus  * J / 2;   // J/2 S+
      sites[i].mpopr(4, m+1, m,   2) =  c_minus * J / 2;   // J/2 S-
    }
    sites[L-1].mpopr(1, m+1, m,   0) =  c_minus * 1.0;     //     S-
    sites[L-1].mpopr(2, m,   m+1, 0) =  c_plus  * 1.0;     //     S+
    mz_plus  -= 1.0;
    mz_minus -= 1.0;
  }
}

void wfn_init(const DmrgInput& input, MpStorages& sites, int nroot)
{
  int L = input.N_sites;
  int M = input.N_max_states;
  int d = input.N_phys_index;
  // set random wfns as initial guess
  if(nroot == 0) {
    sites[0].wfncn.resize(1, d, M);
    for(int i = 1; i < L-1; ++i) {
      sites[i].wfncn.resize(M, d, M);
    }
    sites[L-1].wfncn.resize(M, d, 1);
  }
  for(int i = 0; i < L; ++i) {
    sites[i].wfncn = random_gen;
    Dnormalize(sites[i].wfncn);
  }
}

void str_init(const DmrgInput& input, MpStorages& sites, int nroot)
{
  int L = input.N_sites;
  int M = input.N_max_states;

  sites[L-1].rstrn.resize(1, 1, 1);
  sites[L-1].rstrn = 1.0;
  for(int iroot = 0; iroot < nroot; ++iroot) {
    sites[L-1].rovl0n[iroot].resize(1, 1);
    sites[L-1].rovl0n[iroot] = 1.0;
    sites[L-1].rstr0n[iroot].resize(1, 1, 1);
    sites[L-1].rstr0n[iroot] = 1.0;
  }
  for(int i = L-1; i > 0; --i) {
    canonicalize(0, sites[i].wfncn, sites[i].rmpsn, sites[i-1].wfncn, M);
    sites[i-1].rstrn.free();
    renormalize(0, sites[i].mpopr, sites[i].rstrn, sites[i].rmpsn, sites[i].rmpsn, sites[i-1].rstrn);
    for(int iroot = 0; iroot < nroot; ++iroot) {
      sites[i-1].rovl0n[iroot].free();
      renormalize(0, sites[i].rovl0n[iroot], sites[i].rmpsn, sites[i].rmps0[iroot], sites[i-1].rovl0n[iroot]);
      sites[i-1].rstr0n[iroot].free();
      renormalize(0, sites[i].mpopr, sites[i].rstr0n[iroot], sites[i].rmpsn, sites[i].rmps0[iroot], sites[i-1].rstr0n[iroot]);
    }
  }
  sites[0].lstrn.resize(1, 1, 1);
  sites[0].lstrn = 1.0;
  for(int iroot = 0; iroot < nroot; ++iroot) {
    sites[0].lovl0n[iroot].resize(1, 1);
    sites[0].lovl0n[iroot] = 1.0;
    sites[0].lstr0n[iroot].resize(1, 1, 1);
    sites[0].lstr0n[iroot] = 1.0;
  }
}

double optimize(MpSite& site, const vector<double>& evals, double tole)
{
  int nroot = evals.size();

  vector< DArray<3> > hint0n(nroot, DArray<3>());
  vector< DArray<3> > proj0n(nroot, DArray<3>());
  for(int iroot = 0; iroot < nroot; ++iroot) {
    compute_projector   (site.lovl0n[iroot], site.rovl0n[iroot], site.wfnc0[iroot], proj0n[iroot]);
    compute_sigma_vector(site.mpopr, site.lstr0n[iroot], site.rstr0n[iroot], site.wfnc0[iroot], hint0n[iroot]);
  }
  boost::function<void(const DArray<3>&, DArray<3>&)>
  sgv_functor = boost::bind(compute_sigma_vector, site.mpopr, site.lstrn, site.rstrn, evals, hint0n, proj0n, _1, _2);

  DArray<3> hdiag;
  compute_h_diagonal(site.mpopr, site.lstrn, site.rstrn, hdiag);
  for(int iroot = 0; iroot < nroot; ++iroot) {
    DArray<3>::iterator idiag = hdiag.begin();
    DArray<3>::iterator ihint = hint0n[iroot].begin();
    DArray<3>::iterator iproj = proj0n[iroot].begin();
    double eval = evals[iroot];
    for(; idiag != hdiag.end(); ++idiag, ++ihint, ++iproj) {
      double sdiag = (*iproj) * (eval * (*iproj) - 2.0 * (*ihint));
      *idiag += sdiag;
    }
  }

  double energy = btas::davidson::diagonalize(sgv_functor, hdiag, site.wfncn, tole);

  return energy;
}

double dmrg(ostream& fout, const DmrgInput& input, MpStorages& sites, const vector<double>& evals)
{
  int L = input.N_sites;
  int M = input.N_max_states;
  int d = input.N_phys_index;
  int nroot = evals.size();

  // initialization
  if(nroot == 0) {
    mpo_init(input, sites);
  }
  if(!input.restart) {
    wfn_init(input, sites, nroot);
  }
  str_init(input, sites, nroot);

  int    iter   = 0;
  bool   conv   = false;
  double energy = 0.0;
  double tole   = input.tolerance;

  while(!conv && iter < 100) {
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    fout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    double emin = 1.0e8;
    double wmax = 0.0;
    // fowrad sweep
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fout << "\t\t\tFORWARD SWEEP" << endl;
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(int i = 0; i < L-1; ++i) {
      // diagonalize
      double eiswp = optimize(sites[i], evals, tole / 10);
      if(eiswp < emin) emin = eiswp;
      // canonicalize
      Dcopy(sites[i+1].rmpsn, sites[i+1].wfncn);
      double wiswp = canonicalize(1, sites[i].wfncn, sites[i].lmpsn, sites[i+1].wfncn, M);
      if(wiswp > wmax) wmax = wiswp;
      // print result
      fout.precision(16);
      fout << "\t\t\tEnergy = " << setw(24) << fixed << eiswp << " ( Discarded Weight = ";
      fout.precision(2);
      fout << setw(8) << scientific << wiswp << " ) " << endl;
      // renormalize
      sites[i+1].lstrn.free();
      renormalize(1, sites[i].mpopr, sites[i].lstrn, sites[i].lmpsn, sites[i].lmpsn, sites[i+1].lstrn);
      for(int iroot = 0; iroot < nroot; ++iroot) {
        sites[i+1].lovl0n[iroot].free();
        renormalize(1, sites[i].lovl0n[iroot], sites[i].lmpsn, sites[i].lmps0[iroot], sites[i+1].lovl0n[iroot]);
        sites[i+1].lstr0n[iroot].free();
        renormalize(1, sites[i].mpopr, sites[i].lstr0n[iroot], sites[i].lmpsn, sites[i].lmps0[iroot], sites[i+1].lstr0n[iroot]);
      }
    }
    // backward sweep
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fout << "\t\t\tBACKWARD SWEEP" << endl;
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(int i = L-1; i > 0; --i) {
      // diagonalize
      double eiswp = optimize(sites[i], evals, tole / 10);
      if(eiswp < emin) emin = eiswp;
      // canonicalize
      Dcopy(sites[i-1].lmpsn, sites[i-1].wfncn);
      double wiswp = canonicalize(0, sites[i].wfncn, sites[i].rmpsn, sites[i-1].wfncn, M);
      if(wiswp > wmax) wmax = wiswp;
      // print result
      fout.precision(16);
      fout << "\t\t\tEnergy = " << setw(24) << fixed << eiswp << " ( Discarded Weight = ";
      fout.precision(2);
      fout << setw(8) << scientific << wiswp << " ) " << endl;
      // renormalize
      sites[i-1].rstrn.free();
      renormalize(0, sites[i].mpopr, sites[i].rstrn, sites[i].rmpsn, sites[i].rmpsn, sites[i-1].rstrn);
      for(int iroot = 0; iroot < nroot; ++iroot) {
        sites[i-1].rovl0n[iroot].free();
        renormalize(0, sites[i].rovl0n[iroot], sites[i].rmpsn, sites[i].rmps0[iroot], sites[i-1].rovl0n[iroot]);
        sites[i-1].rstr0n[iroot].free();
        renormalize(0, sites[i].mpopr, sites[i].rstr0n[iroot], sites[i].rmpsn, sites[i].rmps0[iroot], sites[i-1].rstr0n[iroot]);
      }
    }
    // check convergence
    double ediff  = fabs(emin - energy);
    if(ediff < tole) conv = true;
    energy = emin;
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    fout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FISHED" << endl;
    fout.precision(16);
    fout << "\t\t\tEnergy = " << setw(24) << fixed << energy << " ( delta E = ";
    fout.precision(2);
    fout << setw(8) << scientific << ediff << " ) " << endl;
    fout << "\t\t\tLargest discarded weight = " << setw(8) << scientific << wmax << endl;
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    fout << endl;
    ++iter;
  }
  fout.precision(16);
  fout << "\t\t\tFinal Energy = " << setw(24) << fixed << energy   << endl;
  fout << endl;

  return energy;
}

