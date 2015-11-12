#include <iostream>
#include <iomanip>
#include <algorithm>
#include <driver/davidson.h>
#include "driver.h"
#include "driverlrt.h"
#include "dmrg.h"
using namespace std;
using namespace btas;

extern double random_gen(void);

void mpo_init(const DmrgInput& input, MpOprtrs& mpos)
{
  int L = input.N_sites;
  int d = input.N_phys_index;
  mpos.resize(L);
  double J  = input.J;
  double Jz = input.Jz;
  double hz = input.hz;
  double sz = static_cast<double>(d - 1) / 2;
  // resize & set to 0
  mpos[0].resize(1, d, d, 5);
  mpos[0] = 0;
  for(int i = 1; i < L-1; ++i) {
    mpos[i].resize(5, d, d, 5);
    mpos[i] = 0;
  }
  mpos[L-1].resize(5, d, d, 1);
  mpos[L-1] = 0;
  // construct mpos for spin-hamiltonian
  double mz = sz;
  for(int m = 0; m < d; ++m) {
    mpos[0](0, m,   m,   0) =  hz * mz; // hz  Sz
    mpos[0](0, m,   m,   3) =  Jz * mz; // Jz  Sz
    mpos[0](0, m,   m,   4) =  1.0;     // I
    for(int i = 1; i < L-1; ++i) {
      mpos[i](0, m,   m,   0) =  1.0;     //     I
      mpos[i](3, m,   m,   0) =  mz;      //     Sz
      mpos[i](4, m,   m,   0) =  hz * mz; // hz  Sz
      mpos[i](4, m,   m,   3) =  Jz * mz; // Jz  Sz
      mpos[i](4, m,   m,   4) =  1.0;     // I
    }
    mpos[L-1](0, m,   m,   0) =  1.0;     //     I
    mpos[L-1](3, m,   m,   0) =  mz;      //     Sz
    mpos[L-1](4, m,   m,   0) =  hz * mz; // hz  Sz
    mz -= 1.0;
  }
  double mz_plus  = sz - 1.0;
  double mz_minus = sz;
  for(int m = 0; m < d - 1; ++m) {
    double c_plus  = sqrt(sz*(sz+1.0)-mz_plus *(mz_plus +1.0));
    double c_minus = sqrt(sz*(sz+1.0)-mz_minus*(mz_minus-1.0));
    mpos[0](0, m,   m+1, 1) =  c_plus  * J / 2;   // J/2 S+
    mpos[0](0, m+1, m,   2) =  c_minus * J / 2;   // J/2 S-
    for(int i = 1; i < L-1; ++i) {
      mpos[i](1, m+1, m,   0) =  c_minus * 1.0;     //     S-
      mpos[i](2, m,   m+1, 0) =  c_plus  * 1.0;     //     S+
      mpos[i](4, m,   m+1, 1) =  c_plus  * J / 2;   // J/2 S+
      mpos[i](4, m+1, m,   2) =  c_minus * J / 2;   // J/2 S-
    }
    mpos[L-1](1, m+1, m,   0) =  c_minus * 1.0;     //     S-
    mpos[L-1](2, m,   m+1, 0) =  c_plus  * 1.0;     //     S+
    mz_plus  -= 1.0;
    mz_minus -= 1.0;
  }
}

void wfn_init(const DmrgInput& input, MpStates& wfns)
{
  int L = input.N_sites;
  int M = input.N_max_states;
  int d = input.N_phys_index;
  // set random wfns as initial guess
  wfns.resize(L);
  wfns[0].resize(1, d, M);
  wfns[0] = random_gen;
  for(int i = 1; i < L-1; ++i) {
    wfns[i].resize(M, d, M);
    wfns[i] = random_gen;
  }
  wfns[L-1].resize(M, d, 1);
  wfns[L-1] = random_gen;
}

void str_init(const DmrgInput& input,
              const MpOprtrs& mpos, MpStates& wfns,
              DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
              DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr)
{
  int L = input.N_sites;
  int M = input.N_max_states;

  lval.resize(L);
  lmps.resize(L);
  lnul.resize(L);
  lstr.resize(L);

  rval.resize(L);
  rmps.resize(L);
  rnul.resize(L);
  rstr.resize(L);

  rstr[L-1].resize(1, 1, 1);
  rstr[L-1] = 1.0;
  for(int i = L-1; i > 0; --i) {
    canonicalize(0, wfns[i], rmps[i], wfns[i-1], M);
    rstr[i-1].free();
    renormalize (0, mpos[i], rstr[i], rmps[i], rmps[i], rstr[i-1]);
  }

  lstr[ 0 ].resize(1, 1, 1);
  lstr[ 0 ] = 1.0;
}

void str_init(const DmrgInput& input,
              const MpStates& lmps, const vector<MpStates>& lmps0, vector<Overlaps>& lovs,
              const MpStates& rmps, const vector<MpStates>& rmps0, vector<Overlaps>& rovs)
{
  int L = input.N_sites;
  lovs.resize(lmps0.size(), Overlaps(L, DArray<2>()));
  rovs.resize(rmps0.size(), Overlaps(L, DArray<2>()));

  for(int iroot = 0; iroot < rmps0.size(); ++iroot) {
    rovs[iroot][L-1].resize(1, 1);
    rovs[iroot][L-1] = 1.0;
    for(int i = L-1; i > 0; --i) {
      rovs[iroot][i-1].free();
      renormalize (0, rovs[iroot][i], rmps[i], rmps0[iroot][i], rovs[iroot][i-1]);
    }
    lovs[iroot][ 0 ].resize(1, 1);
    lovs[iroot][ 0 ] = 1.0;
  }
}

void dmrg   (ostream& fout, const DmrgInput& input, const MpOprtrs& mpos, MpStates& wfns,
             DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
             DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr)
{
  int L = input.N_sites;
  int M = input.N_max_states;
  int d = input.N_phys_index;

  if(!input.restart) wfn_init(input, wfns);
  str_init(input, mpos, wfns, lval, lmps, lnul, lstr, rval, rmps, rnul, rstr);

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
      // compute diagonal elements
      DArray<3> hdiag;
      compute_h_diagonal(mpos[i], lstr[i], rstr[i], hdiag);
      // set functor for davidson
      boost::function<void(const DArray<3>&, DArray<3>&)>
        sgv_functor = boost::bind(compute_sigma_vector, mpos[i], lstr[i], rstr[i], _1, _2);
      // diagonalize
      double eiswp = btas::davidson::diagonalize(sgv_functor, hdiag, wfns[i], tole / 10);
      if(eiswp < emin) emin = eiswp;
      // compute density matrix ( l, n, l', n' )
      DArray<4> lsdm;
      Dgemm(NoTrans, Trans, 1.0, wfns[i], wfns[i], 1.0, lsdm);
      // canonicalize
      double wiswp = canonicalize(1, lsdm, lval[i], lmps[i], lnul[i], M);
      if(wiswp > wmax) wmax = wiswp;
      // compute guess wfn for next
      DArray<2> resm;
      Dgemm(Trans, NoTrans, 1.0, lmps[i], wfns[i], 1.0, resm);
      wfns[i+1].free();
      Dgemm(NoTrans, NoTrans, 1.0, resm, rmps[i+1], 1.0, wfns[i+1]);
      // print result
      fout.precision(16);
      fout << "\t\t\tEnergy = " << setw(24) << fixed << eiswp << " ( Discarded Weight = ";
      fout.precision(2);
      fout << setw(8) << scientific << wiswp << " ) " << endl;
      // renormalize
      lstr[i+1].free();
      renormalize(1, mpos[i], lstr[i], lmps[i], lmps[i], lstr[i+1]);
    }
    // backward sweep
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fout << "\t\t\tBACKWARD SWEEP" << endl;
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(int i = L-1; i > 0; --i) {
      // compute diagonal elements
      DArray<3> hdiag;
      compute_h_diagonal(mpos[i], lstr[i], rstr[i], hdiag);
      // set functor for davidson
      boost::function<void(const DArray<3>&, DArray<3>&)>
        sgv_functor = boost::bind(compute_sigma_vector, mpos[i], lstr[i], rstr[i], _1, _2);
      // diagonalize
      double eiswp = btas::davidson::diagonalize(sgv_functor, hdiag, wfns[i], tole / 10);
      if(eiswp < emin) emin = eiswp;
      // compute density matrix ( n, r, n', r' )
      DArray<4> rsdm;
      Dgemm(Trans, NoTrans, 1.0, wfns[i], wfns[i], 1.0, rsdm);
      // canonicalize
      double wiswp = canonicalize(0, rsdm, rval[i], rmps[i], rnul[i], M);
      if(wiswp > wmax) wmax = wiswp;
      // compute guess wfn for next
      DArray<2> resm;
      Dgemm(NoTrans, Trans, 1.0, wfns[i], rmps[i], 1.0, resm);
      wfns[i-1].free();
      Dgemm(NoTrans, NoTrans, 1.0, lmps[i-1], resm, 1.0, wfns[i-1]);
      // print result
      fout.precision(16);
      fout << "\t\t\tEnergy = " << setw(24) << fixed << eiswp << " ( Discarded Weight = ";
      fout.precision(2);
      fout << setw(8) << scientific << wiswp << " ) " << endl;
      // renormalize
      rstr[i-1].free();
      renormalize(0, mpos[i], rstr[i], rmps[i], rmps[i], rstr[i-1]);
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
}

void dmrg   (ostream& fout, const DmrgInput& input, const MpOprtrs& mpos, MpStates& wfns,
             DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
             DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr,
             const vector<MpStates>& wfns0, const vector<MpStates>& lmps0, const vector<MpStates>& rmps0)
{
  int L = input.N_sites;
  int M = input.N_max_states;
  int d = input.N_phys_index;

  if(!input.restart) {
    wfns.resize(L);
    for(int i = 0; i < L; ++i) {
      wfns[i].resize(wfns0[0][i].shape());
      wfns[i] = random_gen;
    }
  }
  str_init(input, mpos, wfns, lval, lmps, lnul, lstr, rval, rmps, rnul, rstr);

  vector< Overlaps > lovs;
  vector< Overlaps > rovs;
  str_init(input, lmps, lmps0, lovs, rmps, rmps0, rovs);

  int    nroot  = wfns0.size();

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
      // compute diagonal elements
      DArray<3> hdiag;
      compute_h_diagonal(mpos[i], lstr[i], rstr[i], hdiag);
      // set functor for davidson
      boost::function<void(const DArray<3>&, DArray<3>&)>
        sgv_functor = boost::bind(compute_sigma_vector, mpos[i], lstr[i], rstr[i], _1, _2);
      // compute overlaps
      vector< DArray<3> > ovlpn(nroot, DArray<3>());
      vector< DArray<3> > wfn0n(nroot, DArray<3>());
      for(int iroot = 0; iroot < nroot; ++iroot) {
        DArray<3> scr1;
        Dgemm(NoTrans, NoTrans, 1.0, lovs[iroot][i], wfns0[iroot][i], 1.0, scr1);
        Dgemm(NoTrans,   Trans, 1.0, scr1, rovs[iroot][i], 1.0, ovlpn[iroot]);
        Dcopy(wfns0[iroot][i], wfn0n[iroot]);
      }
      // diagonalize
      double eiswp = btas::davidson::diagonalize(sgv_functor, hdiag, wfns[i], ovlpn, wfn0n, tole / 10);
// debug
      for(int iroot = 0; iroot < nroot; ++iroot) {
        double iovlp = Ddot(wfns[i], ovlpn[iroot]);
        cout << "\t\toverlap[" << setw(2) << iroot << "] = " << iovlp << endl;
      }
// debug
      if(eiswp < emin) emin = eiswp;
      // compute density matrix ( l, n, l', n' )
      DArray<4> lsdm;
      Dgemm(NoTrans, Trans, 1.0, wfns[i], wfns[i], 1.0, lsdm);
      // canonicalize
      double wiswp = canonicalize(1, lsdm, lval[i], lmps[i], lnul[i], M);
      if(wiswp > wmax) wmax = wiswp;
      // compute guess wfn for next
      DArray<2> resm;
      Dgemm(Trans, NoTrans, 1.0, lmps[i], wfns[i], 1.0, resm);
      wfns[i+1].free();
      Dgemm(NoTrans, NoTrans, 1.0, resm, rmps[i+1], 1.0, wfns[i+1]);
      // print result
      fout.precision(16);
      fout << "\t\t\tEnergy = " << setw(24) << fixed << eiswp << " ( Discarded Weight = ";
      fout.precision(2);
      fout << setw(8) << scientific << wiswp << " ) " << endl;
      // renormalize
      lstr[i+1].free();
      renormalize(1, mpos[i], lstr[i], lmps[i], lmps[i], lstr[i+1]);
      for(int iroot = 0; iroot < nroot; ++iroot) {
        lovs[iroot][i+1].free();
        renormalize(1, lovs[iroot][i], lmps[i], lmps0[iroot][i], lovs[iroot][i+1]);
      }
    }
    // backward sweep
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fout << "\t\t\tBACKWARD SWEEP" << endl;
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(int i = L-1; i > 0; --i) {
      // compute diagonal elements
      DArray<3> hdiag;
      compute_h_diagonal(mpos[i], lstr[i], rstr[i], hdiag);
      // set functor for davidson
      boost::function<void(const DArray<3>&, DArray<3>&)>
        sgv_functor = boost::bind(compute_sigma_vector, mpos[i], lstr[i], rstr[i], _1, _2);
      // compute overlaps
      vector< DArray<3> > ovlpn(nroot, DArray<3>());
      vector< DArray<3> > wfn0n(nroot, DArray<3>());
      for(int iroot = 0; iroot < nroot; ++iroot) {
        DArray<3> scr1;
        Dgemm(NoTrans, NoTrans, 1.0, lovs[iroot][i], wfns0[iroot][i], 1.0, scr1);
        Dgemm(NoTrans,   Trans, 1.0, scr1, rovs[iroot][i], 1.0, ovlpn[iroot]);
        Dcopy(wfns0[iroot][i], wfn0n[iroot]);
      }
      // diagonalize
      double eiswp = btas::davidson::diagonalize(sgv_functor, hdiag, wfns[i], ovlpn, wfn0n, tole / 10);
// debug
      for(int iroot = 0; iroot < nroot; ++iroot) {
        double iovlp = Ddot(wfns[i], ovlpn[iroot]);
        cout << "\t\toverlap[" << setw(2) << iroot << "] = " << iovlp << endl;
      }
// debug
      if(eiswp < emin) emin = eiswp;
      // compute density matrix ( n, r, n', r' )
      DArray<4> rsdm;
      Dgemm(Trans, NoTrans, 1.0, wfns[i], wfns[i], 1.0, rsdm);
      // canonicalize
      double wiswp = canonicalize(0, rsdm, rval[i], rmps[i], rnul[i], M);
      if(wiswp > wmax) wmax = wiswp;
      // compute guess wfn for next
      DArray<2> resm;
      Dgemm(NoTrans, Trans, 1.0, wfns[i], rmps[i], 1.0, resm);
      wfns[i-1].free();
      Dgemm(NoTrans, NoTrans, 1.0, lmps[i-1], resm, 1.0, wfns[i-1]);
      // print result
      fout.precision(16);
      fout << "\t\t\tEnergy = " << setw(24) << fixed << eiswp << " ( Discarded Weight = ";
      fout.precision(2);
      fout << setw(8) << scientific << wiswp << " ) " << endl;
      // renormalize
      rstr[i-1].free();
      renormalize(0, mpos[i], rstr[i], rmps[i], rmps[i], rstr[i-1]);
      for(int iroot = 0; iroot < nroot; ++iroot) {
        rovs[iroot][i-1].free();
        renormalize(0, rovs[iroot][i], rmps[i], rmps0[iroot][i], rovs[iroot][i-1]);
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
}

void dmrglrt(ostream& fout, const DmrgInput& input, const MpOprtrs& mpos, MpStates& wfns,
             DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
             DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr, int nroot)
{
  int L = input.N_sites;
  int M = input.N_max_states;
  int d = input.N_phys_index;

  // allocate working spaces and initialize
  vector< vector< DArray<3> > > wfns1st(L, vector< DArray<3> >(nroot, DArray<3>()));
  vector< vector< DArray<3> > > lmps1st(L, vector< DArray<3> >(nroot, DArray<3>()));
  vector< vector< DArray<3> > > rmps1st(L, vector< DArray<3> >(nroot, DArray<3>()));
  vector< vector< DArray<3> > > lstr1st(L, vector< DArray<3> >(nroot, DArray<3>()));
  vector< vector< DArray<3> > > rstr1st(L, vector< DArray<3> >(nroot, DArray<3>()));
  for(int i = 0; i < L; ++i) {
    for(int iroot = 0; iroot < nroot; ++iroot) {
      wfns1st[i][iroot].resize(wfns[i].shape()); wfns1st[i][iroot] = random_gen;
      lmps1st[i][iroot].resize(lmps[i].shape()); lmps1st[i][iroot] = 0.0;
      rmps1st[i][iroot].resize(rmps[i].shape()); rmps1st[i][iroot] = 0.0;
      lstr1st[i][iroot].resize(lstr[i].shape()); lstr1st[i][iroot] = 0.0;
      rstr1st[i][iroot].resize(rstr[i].shape()); rstr1st[i][iroot] = 0.0;
    }
  }

  int    iter   = 0;
  bool   conv   = false;
  double tole   = input.tolerance;

  vector<double> energy(nroot, 0.0);

  while(!conv && iter < 100) {
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    fout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    vector<double> emin(nroot, 1.0e8);
    // fowrad sweep
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fout << "\t\t\tFORWARD SWEEP" << endl;
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(int i = 0; i < L-1; ++i) {
      // compute diagonal elements
      DArray<3> hdiag;
      compute_h_diagonal(mpos[i], lstr[i], rstr[i], hdiag);
      // diagonalize
      vector<double> eiswp = lrt::diagonalize(mpos[i], hdiag, lstr[i], rstr[i], wfns[i],
                                              lstr1st[i], rstr1st[i], wfns1st[i], nroot, tole / 10);
      for(int iroot = 0; iroot < nroot; ++iroot) {
        if(eiswp[iroot] < emin[iroot]) emin[iroot] = eiswp[iroot];
        // compute 1st-order density matrix ( l, n, l', n' )
        DArray<4> lsdm1st;
        Dgemm(NoTrans, Trans, 1.0, wfns1st[i][iroot], wfns[i], 1.0, lsdm1st);
        Dgemm(NoTrans, Trans, 1.0, wfns[i], wfns1st[i][iroot], 1.0, lsdm1st);
        // canonicalize to get 1st-order rotation matrix
        lrt::canonicalize(1, lsdm1st, lval[i], lmps[i], lnul[i], lmps1st[i][iroot]);
        // print result
        fout.precision(16);
        fout << "\t\t\tEnergy[" << setw(2) << iroot << "] = " << setw(24) << fixed << eiswp[iroot] << endl;
        // renormalize
        lstr1st[i+1][iroot].free();
        renormalize(1, mpos[i], lstr1st[i][iroot], lmps[i], lmps[i], lstr1st[i+1][iroot]);
        renormalize(1, mpos[i], lstr[i], lmps1st[i][iroot], lmps[i], lstr1st[i+1][iroot]);
        renormalize(1, mpos[i], lstr[i], lmps[i], lmps1st[i][iroot], lstr1st[i+1][iroot]);
      }
      fout << "\t--------------------------------------------------"
           <<   "--------------------------------------------------" << endl;
    }
    // backward sweep
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fout << "\t\t\tBACKWARD SWEEP" << endl;
    fout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++"
         <<   "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(int i = L-1; i > 0; --i) {
      // compute diagonal elements
      DArray<3> hdiag;
      compute_h_diagonal(mpos[i], lstr[i], rstr[i], hdiag);
      // diagonalize
      vector<double> eiswp = lrt::diagonalize(mpos[i], hdiag, lstr[i], rstr[i], wfns[i],
                                              lstr1st[i], rstr1st[i], wfns1st[i], nroot, tole / 10);
      for(int iroot = 0; iroot < nroot; ++iroot) {
        if(eiswp[iroot] < emin[iroot]) emin[iroot] = eiswp[iroot];
        // compute 1st-order density matrix ( l, n, l', n' )
        DArray<4> rsdm1st;
        Dgemm(Trans, NoTrans, 1.0, wfns1st[i][iroot], wfns[i], 1.0, rsdm1st);
        Dgemm(Trans, NoTrans, 1.0, wfns[i], wfns1st[i][iroot], 1.0, rsdm1st);
        // canonicalize to get 1st-order rotation matrix
        lrt::canonicalize(0, rsdm1st, rval[i], rmps[i], rnul[i], rmps1st[i][iroot]);
        // print result
        fout.precision(16);
        fout << "\t\t\tEnergy[" << setw(2) << iroot << "] = " << setw(24) << fixed << eiswp[iroot] << endl;
        // renormalize
        rstr1st[i-1][iroot].free();
        renormalize(0, mpos[i], rstr1st[i][iroot], rmps[i], rmps[i], rstr1st[i-1][iroot]);
        renormalize(0, mpos[i], rstr[i], rmps1st[i][iroot], rmps[i], rstr1st[i-1][iroot]);
        renormalize(0, mpos[i], rstr[i], rmps[i], rmps1st[i][iroot], rstr1st[i-1][iroot]);
      }
      fout << "\t--------------------------------------------------"
           <<   "--------------------------------------------------" << endl;
    }
    // check convergence
    conv = true;
    for(int iroot = 0; iroot < nroot; ++iroot) {
      double ediff  = fabs(emin[iroot] - energy[iroot]);
      if(ediff >= tole) conv = false;
    }
    energy = emin;
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    fout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FISHED" << endl;
    fout.precision(16);
    for(int iroot = 0; iroot < nroot; ++iroot) {
      fout << "\t\t\tEnergy[" << setw(2) << iroot << "] = " << setw(24) << fixed << energy[iroot] << endl;
    }
    fout << "\t=================================================="
         <<   "==================================================" << endl;
    fout << endl;
    ++iter;
  }
}
