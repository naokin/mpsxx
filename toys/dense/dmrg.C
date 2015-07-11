#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <legacy/DArray.h>
#include <legacy/Dblas.h>
#include "btas_template_specialize.h"

#include "dmrg.h"
#include "lrt.h"
#include "driver.h"
#include "davidson.h"
using namespace btas;

//
// random number generator
//
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; }

//
// Heisenberg model
//
void prototype::Heisenberg::construct_mpo(MpStorages& sites, const HeisenbergModel& params)
{
  int    Nz = params.Nz;
  double J  = params.J;
  double Jz = params.Jz;
  double Hz = params.Hz;

  int    L  = sites.size();
  int    d  = Nz + 1;
  double sz = static_cast<double>(Nz)/2;

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tCONSTRUCT MATRIX PRODUCT OPERATORS (MPOs) "                                                         << endl;
  cout.precision(4);
  cout << "\t\t\t+ coupling coefficient J  : " << setw(8) << fixed << J  << endl; 
  cout << "\t\t\t+ coupling coefficient Jz : " << setw(8) << fixed << Jz << endl;
  cout << "\t\t\t+ coupling coefficient Hz : " << setw(8) << fixed << Hz << endl;
  cout << "\t\t\t+ coupling coefficient Sz : " << setw(8) << fixed << sz << endl;
  cout << "\t====================================================================================================" << endl;

  // resize & set to 0
  sites[ 0 ].mpo.resize(1, d, d, 5);
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpo.resize(5, d, d, 5);
  }
  sites[L-1].mpo.resize(5, d, d, 1);

  // construct mpos for spin-hamiltonian
  double mz = sz;
  for(int m = 0; m < d; ++m) {
    double data_Id = 1.0;
    double data_Mz = mz;
    double data_Hz = Hz * mz;
    double data_Jz = Jz * mz;

    sites[ 0 ].mpo(0, m, m, 0) = data_Hz; // Hz  Sz
    sites[ 0 ].mpo(0, m, m, 3) = data_Jz; // Jz  Sz
    sites[ 0 ].mpo(0, m, m, 4) = data_Id; //     I
    for(int i = 1; i < L-1; ++i) {
      sites[i].mpo(0, m, m, 0) = data_Id; //     I
      sites[i].mpo(3, m, m, 0) = data_Mz; //     Sz
      sites[i].mpo(4, m, m, 0) = data_Hz; // Hz  Sz
      sites[i].mpo(4, m, m, 3) = data_Jz; // Jz  Sz
      sites[i].mpo(4, m, m, 4) = data_Id; //     I
    }
    sites[L-1].mpo(0, m, m, 0) = data_Id; //     I
    sites[L-1].mpo(3, m, m, 0) = data_Mz; //     Sz
    sites[L-1].mpo(4, m, m, 0) = data_Hz; // Hz  Sz
    mz -= 1.0;
  }
  double mz_plus  = sz - 1.0;
  double mz_minus = sz;
  for(int m = 0; m < d - 1; ++m) {
    double c_plus  = sqrt(sz*(sz+1.0)-mz_plus *(mz_plus +1.0));
    double c_minus = sqrt(sz*(sz+1.0)-mz_minus*(mz_minus-1.0));
    // set block elements
    double data_Sp = c_plus;
    double data_Sm = c_minus;
    double data_Jp = c_plus  * J / 2;
    double data_Jm = c_minus * J / 2;
    // insert blocks
    sites[ 0 ].mpo(0, m,   m+1, 1) = data_Jp; // J/2 S+
    sites[ 0 ].mpo(0, m+1, m,   2) = data_Jm; // J/2 S-
    for(int i = 1; i < L-1; ++i) {
      sites[i].mpo(1, m+1, m,   0) = data_Sm; //     S-
      sites[i].mpo(2, m,   m+1, 0) = data_Sp; //     S+
      sites[i].mpo(4, m,   m+1, 1) = data_Jp; // J/2 S+
      sites[i].mpo(4, m+1, m,   2) = data_Jm; // J/2 S-
    }
    sites[L-1].mpo(1, m+1, m,   0) = data_Sm; //     S-
    sites[L-1].mpo(2, m,   m+1, 0) = data_Sp; //     S+
    mz_plus  -= 1.0;
    mz_minus -= 1.0;
  }
}

//
// Hubbard model
//
void prototype::Hubbard::construct_mpo(MpStorages& sites, const HubbardModel& params)
{
  double t  = params.t;
  double U  = params.U;

  int    L  = sites.size();

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tCONSTRUCT MATRIX PRODUCT OPERATORS (MPOs) "                                                         << endl;
  cout.precision(4);
  cout << "\t\t\t+ coupling coefficient t  : " << setw(8) << fixed << t  << endl; 
  cout << "\t\t\t+ coupling coefficient U  : " << setw(8) << fixed << U  << endl;
  cout << "\t====================================================================================================" << endl;

  // resize & set to 0
  sites[ 0 ].mpo.resize(1, 4, 4, 6);
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpo.resize(6, 4, 4, 6);
  }
  sites[L-1].mpo.resize(6, 4, 4, 1);

  // set block elements
  double data_Ip = 1.0;
  double data_Im =-1.0;
  double data_tp = t;
  double data_tm =-t;
  double data_Un = U;
  // insert blocks
  sites[ 0 ].mpo(0, 3, 3, 0) = data_Un; //  U
  sites[ 0 ].mpo(0, 1, 0, 1) = data_tp; //  t a+
  sites[ 0 ].mpo(0, 3, 2, 1) = data_tm; //  t a+ [x(-1) due to P(a,b)]
  sites[ 0 ].mpo(0, 0, 1, 2) = data_tm; // -t a-
  sites[ 0 ].mpo(0, 2, 3, 2) = data_tp; // -t a- [x(-1) due to P(a,b)]
  sites[ 0 ].mpo(0, 2, 0, 3) = data_tp; //  t b+
  sites[ 0 ].mpo(0, 3, 1, 3) = data_tp; //  t b+
  sites[ 0 ].mpo(0, 0, 2, 4) = data_tm; // -t b-
  sites[ 0 ].mpo(0, 1, 3, 4) = data_tm; // -t b-
  sites[ 0 ].mpo(0, 0, 0, 5) = data_Ip; //  I
  sites[ 0 ].mpo(0, 1, 1, 5) = data_Ip; //  I
  sites[ 0 ].mpo(0, 2, 2, 5) = data_Ip; //  I
  sites[ 0 ].mpo(0, 3, 3, 5) = data_Ip; //  I
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpo(0, 0, 0, 0) = data_Ip; //  I
    sites[i].mpo(0, 1, 1, 0) = data_Ip; //  I
    sites[i].mpo(0, 2, 2, 0) = data_Ip; //  I
    sites[i].mpo(0, 3, 3, 0) = data_Ip; //  I
    sites[i].mpo(1, 0, 1, 0) = data_Ip; //  a-
    sites[i].mpo(1, 2, 3, 0) =-data_Im; //  a- [x(-1) due to P(a,b)]
    sites[i].mpo(2, 1, 0, 0) =-data_Ip; //  a+
    sites[i].mpo(2, 3, 2, 0) = data_Im; //  a+ [x(-1) due to P(a,b)]
    sites[i].mpo(3, 0, 2, 0) = data_Ip; //  b-
    sites[i].mpo(3, 1, 3, 0) =-data_Ip; //  b-
    sites[i].mpo(4, 2, 0, 0) =-data_Ip; //  b+
    sites[i].mpo(4, 3, 1, 0) = data_Ip; //  b+
    sites[i].mpo(5, 3, 3, 0) = data_Un; //  U
    sites[i].mpo(5, 1, 0, 1) = data_tp; //  t a+
    sites[i].mpo(5, 3, 2, 1) = data_tm; //  t a+ [x(-1) due to P(a,b)]
    sites[i].mpo(5, 0, 1, 2) = data_tm; // -t a-
    sites[i].mpo(5, 2, 3, 2) = data_tp; // -t a- [x(-1) due to P(a,b)]
    sites[i].mpo(5, 2, 0, 3) = data_tp; //  t b+
    sites[i].mpo(5, 3, 1, 3) = data_tp; //  t b+
    sites[i].mpo(5, 0, 2, 4) = data_tm; // -t b-
    sites[i].mpo(5, 1, 3, 4) = data_tm; // -t b-
    sites[i].mpo(5, 0, 0, 5) = data_Ip; //  I
    sites[i].mpo(5, 1, 1, 5) = data_Ip; //  I
    sites[i].mpo(5, 2, 2, 5) = data_Ip; //  I
    sites[i].mpo(5, 3, 3, 5) = data_Ip; //  I
  }
  sites[L-1].mpo(0, 0, 0, 0) = data_Ip; //  I
  sites[L-1].mpo(0, 1, 1, 0) = data_Ip; //  I
  sites[L-1].mpo(0, 2, 2, 0) = data_Ip; //  I
  sites[L-1].mpo(0, 3, 3, 0) = data_Ip; //  I
  sites[L-1].mpo(1, 0, 1, 0) = data_Ip; //  a-
  sites[L-1].mpo(1, 2, 3, 0) =-data_Im; //  a- [x(-1) due to P(a,b)]
  sites[L-1].mpo(2, 1, 0, 0) =-data_Ip; //  a+
  sites[L-1].mpo(2, 3, 2, 0) = data_Im; //  a+ [x(-1) due to P(a,b)]
  sites[L-1].mpo(3, 0, 2, 0) = data_Ip; //  b-
  sites[L-1].mpo(3, 1, 3, 0) =-data_Ip; //  b-
  sites[L-1].mpo(4, 2, 0, 0) =-data_Ip; //  b+
  sites[L-1].mpo(4, 3, 1, 0) = data_Ip; //  b+
  sites[L-1].mpo(5, 3, 3, 0) = data_Un; //  U
}

void prototype::initialize(MpStorages& sites, int M)
{
  int L = sites.size();

  int Mx = 1;
  int d;

  for(int i = 0; i < L-1; ++i) {
    d = sites[i].mpo.extent(2);
    sites[i].wfnc[0].resize(Mx, d, M);
    sites[i].wfnc[0] = rgen;
    sites[i].save(i);
    Mx = M;
  }

  d = sites[L-1].mpo.extent(2);
  sites[L-1].wfnc[0].resize(Mx, d, 1);
  sites[L-1].wfnc[0] = rgen;

  //
  // canonicalize & renormalize
  //

  sites[L-1].ropr[0].resize(1, 1, 1);
  sites[L-1].ropr[0] = 1.0;

  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    util::Normalize(sites[i].wfnc[0]);
    Canonicalize(0, sites[i].wfnc[0], sites[i].rmps[0], M);
    Dcopy(sites[i-1].wfnc[0], sites[i-1].lmps[0]);
    sites[i-1].wfnc[0].free();
    ComputeGuess(0, sites[i].rmps[0], sites[i].wfnc[0], sites[i-1].lmps[0], sites[i-1].wfnc[0]);
    sites[i-1].ropr[0].free();
    Renormalize (0, sites[i].mpo, sites[i].ropr[0], sites[i].rmps[0], sites[i].rmps[0], sites[i-1].ropr[0]);
    sites[i].save(i);
  }

  util::Normalize(sites[0].wfnc[0]);
  sites[0].lopr[0].resize(1, 1, 1);
  sites[0].lopr[0] = 1.0;
  sites[0].save(0);
}

double prototype::optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M, double noise)
{
  boost::function<void(const DArray<3>&, DArray<3>&)>
  f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, sysdot.lopr[0], sysdot.ropr[0], _1, _2);
  DArray<3> diag;
  ComputeDiagonal(sysdot.mpo, sysdot.lopr[0], sysdot.ropr[0], diag);
  double energy = davidson::diagonalize(f_contract, diag, sysdot.wfnc[0]);

  if(forward) {
    Canonicalize(1, sysdot.wfnc[0], sysdot.lmps[0], sysdot.ltnj[0]);
//  Canonicalize(1, sysdot.wfnc[0], sysdot.lmps[0], M);
    envdot.wfnc[0].free();
    ComputeGuess(1, sysdot.lmps[0], sysdot.wfnc[0], envdot.rmps[0], envdot.wfnc[0]);
    envdot.lopr[0].free();
    Renormalize (1, sysdot.mpo,  sysdot.lopr[0], sysdot.lmps[0], sysdot.lmps[0], envdot.lopr[0]);
  }
  else {
    Canonicalize(0, sysdot.wfnc[0], sysdot.rmps[0], sysdot.rtnj[0]);
//  Canonicalize(0, sysdot.wfnc[0], sysdot.rmps[0], M);
    envdot.wfnc[0].free();
    ComputeGuess(0, sysdot.rmps[0], sysdot.wfnc[0], envdot.lmps[0], envdot.wfnc[0]);
    envdot.ropr[0].free();
    Renormalize (0, sysdot.mpo,  sysdot.ropr[0], sysdot.rmps[0], sysdot.rmps[0], envdot.ropr[0]);
  }

  return energy;
}

double prototype::optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M, double noise)
{
  DArray<4> wfnc;
  DArray<4> diag;
  boost::function<void(const DArray<4>&, DArray<4>&)> f_contract;
  if(forward) {
    Dgemm(NoTrans, NoTrans, 1.0, sysdot.wfnc[0], envdot.rmps[0], 1.0, wfnc);
    f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, envdot.mpo, sysdot.lopr[0], envdot.ropr[0], _1, _2);
    ComputeDiagonal(sysdot.mpo, envdot.mpo, sysdot.lopr[0], envdot.ropr[0], diag);
  }
  else {
    Dgemm(NoTrans, NoTrans, 1.0, envdot.lmps[0], sysdot.wfnc[0], 1.0, wfnc);
    f_contract = boost::bind(ComputeSigmaVector, envdot.mpo, sysdot.mpo, envdot.lopr[0], sysdot.ropr[0], _1, _2);
    ComputeDiagonal(envdot.mpo, sysdot.mpo, envdot.lopr[0], sysdot.ropr[0], diag);
  }

  double energy = davidson::diagonalize(f_contract, diag, wfnc);

  if(forward) {
    Canonicalize(1,        wfnc, sysdot.lmps[0], envdot.wfnc[0], M);
    envdot.lopr[0].free();
    Renormalize (1, sysdot.mpo,  sysdot.lopr[0], sysdot.lmps[0], sysdot.lmps[0], envdot.lopr[0]);
  }
  else {
    Canonicalize(0,        wfnc, sysdot.rmps[0], envdot.wfnc[0], M);
    envdot.ropr[0].free();
    Renormalize (0, sysdot.mpo,  sysdot.ropr[0], sysdot.rmps[0], sysdot.rmps[0], envdot.ropr[0]);
  }

  return energy;
}

double prototype::dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M, double noise)
{
  int    L    = sites.size();
  double emin = 1.0e8;
  // fowrad sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  sites[0].load(0);
  for(int i = 0; i < L-1; ++i) {
    sites[i+1].load(i+1);
    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(1, sites[i], sites[i+1], M, noise);
    else                eswp = optimize_twosite(1, sites[i], sites[i+1], M, noise);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    sites[i].save(i);
  }
  // backward sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(0, sites[i], sites[i-1], M, noise);
    else                eswp = optimize_twosite(0, sites[i], sites[i-1], M, noise);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    sites[i].save(i);
  }
  sites[0].save(0);
  return emin;
}

double prototype::dmrg(const DmrgInput& input, MpStorages& sites)
{
  const int L = input.N_sites;
  const int M = input.N_max_states;

  sites.resize(L, MpSite("state-0", 1, input.prefix));

  switch(input.model) {
    case HEISENBERG:
      Heisenberg::construct_mpo(sites, input.heisenberg);
      initialize(sites, M);
      break;
    case HUBBARD:
      Hubbard::construct_mpo(sites, input.hubbard);
      initialize(sites, M);
      break;
    default:
      cout << "\tUnknown model type specified" << endl;
      return 0.0;
  }

  double esav = 1.0e8;

  //
  // two-site optimization
  //
  for(int iter = 0; iter < 100; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;

    double noise;
    if(iter < 10)
      noise = 1.0e-4;
    else if(iter < 20)
      noise = 1.0e-6;
    else
      noise = 0.0;

    double eswp = dmrg_sweep(sites, TWOSITE, M, noise);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(fabs(edif) < 1.0e-8) break;
  }

  //
  // one-site optimization
  //
  for(int iter = 0; iter < 100; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(sites, ONESITE, M, 0.0);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(fabs(edif) < 1.0e-8) break;
  }

  construct_effective_hamiltonian(sites, esav);

  return esav;
}
