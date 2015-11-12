#include <iostream>
#include <iomanip>
#include <cmath>

#include <legacy/QSPARSE/QSTArray.h>

#include "make_sweep.h"
#include "fermion.h"
#include "fileio.h"
#include "optimize.hpp"

double mpsxx::make_sweep(const std::vector<double>& E, const mpsxx::DMRGInput& input, const size_t& iroot)
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  size_t N = input.N_sites;
  int    M = input.N_max_states;
  double T = input.tolerance;

  btas::QSTArray<double,4,fermion> OpSys;
  btas::QSTArray<double,4,fermion> OpEnv;

  std::vector<btas::QSTArray<double,3,fermion>> WfSys(iroot+1);
  std::vector<btas::QSTArray<double,3,fermion>> WfEnv(iroot+1);
  std::vector<btas::QSTArray<double,3,fermion>> WfMps(iroot+1);

  std::vector<btas::QSTArray<double,3,fermion>> lHopr(iroot+1);
  std::vector<btas::QSTArray<double,3,fermion>> rHopr(iroot+1);

  std::vector<btas::QSTArray<double,2,fermion>> lSopr(iroot+1);
  std::vector<btas::QSTArray<double,2,fermion>> rSopr(iroot+1);

  double emin = 1.0e8;

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  load(OpSys,getfile("mpo",input.prefix,0));
  for(size_t k = 0; k <= iroot; ++k) {
    load(WfSys[k],getfile("wave",  input.prefix,0,k));
    load(lHopr[k],getfile("left-H",input.prefix,0,iroot,k));
    load(lSopr[k],getfile("left-S",input.prefix,0,iroot,k));
  }

  for(size_t i = 0; i < N-1; ++i) {

    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSITE [ " << setw(3) << i << " ] " << endl;
    cout << "\t----------------------------------------------------------------------------------------------------" << endl;

    cout << "\t\tloading operators and wavefunction of next site (env)" << endl;
    load(OpEnv,getfile("mpo",input.prefix,i+1));
    for(size_t k = 0; k <= iroot; ++k) {
      load(WfEnv[k],getfile("rmps",input.prefix,i+1,k));
    }
    for(size_t k = 0; k < iroot; ++k) {
      load(WfMps[k],getfile("lmps",input.prefix,i,k));
    }

    double eswp;

    // Optimize for each site
    if(input.algorithm == ONESITE) {
      cout << "\t\toptimizing wavefunction: 1-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(rHopr[k],getfile("right-H",input.prefix,i,iroot,k));
        load(rSopr[k],getfile("right-S",input.prefix,i,iroot,k));
      }
//    eswp = optimize_onesite_merged(1,E,OpSys,      lHopr,rHopr,lSopr,rSopr,WfSys,WfEnv,0.1*T,M);
      eswp = optimize_onesite(1,E,OpSys,      lHopr,rHopr,lSopr,rSopr,WfMps,WfSys,WfEnv,0.1*T,M);
    }
    else {
      cout << "\t\toptimizing wavefunction: 2-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(rHopr[k],getfile("right-H",input.prefix,i+1,iroot,k));
        load(rSopr[k],getfile("right-S",input.prefix,i+1,iroot,k));
      }
//    eswp = optimize_twosite_merged(1,E,OpSys,OpEnv,lHopr,rHopr,lSopr,rSopr,WfSys,WfEnv,0.1*T,M);
      eswp = optimize_twosite(1,E,OpSys,OpEnv,lHopr,rHopr,lSopr,rSopr,WfMps,WfSys,WfEnv,0.1*T,M);
    }
    if(eswp < emin) emin = eswp;
    // print result
    cout << "\t\t--------------------------------------------------------------------------------" << endl;
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    cout << "\t\t--------------------------------------------------------------------------------" << endl;

    cout << "\t\tsaving operators and wavefunction of this site (sys)" << endl;
    save(WfSys[iroot],getfile("lmps",input.prefix,i,iroot));
    for(size_t k = 0; k <= iroot; ++k) {
      save(lHopr[k],getfile("left-H",input.prefix,i+1,iroot,k));
      save(lSopr[k],getfile("left-S",input.prefix,i+1,iroot,k));
    }
    // load WfSys for the next loop
    for(size_t k = 0; k < iroot; ++k) {
      load(WfSys[k],getfile("wave",input.prefix,i+1,k));
    }
    OpSys = OpEnv;
    WfSys[iroot] = WfEnv[iroot];
    save(WfSys[iroot],getfile("wave",input.prefix,i+1,iroot));
  }

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  for(size_t k = 0; k <= iroot; ++k) {
    load(rHopr[k],getfile("right-H",input.prefix,N-1,iroot,k));
    load(rSopr[k],getfile("right-S",input.prefix,N-1,iroot,k));
  }

  for(size_t i = N-1; i > 0; --i) {

    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSITE [ " << setw(3) << i << " ] " << endl;
    cout << "\t----------------------------------------------------------------------------------------------------" << endl;

    cout << "\t\tloading operators and wavefunction of next site (env)" << endl;
    load(OpEnv,getfile("mpo",input.prefix,i-1));
    for(size_t k = 0; k <= iroot; ++k) {
      load(WfEnv[k],getfile("lmps",input.prefix,i-1,k));
    }
    for(size_t k = 0; k < iroot; ++k) {
      load(WfMps[k],getfile("rmps",input.prefix,i,k));
    }

    double eswp;

    // Optimize for each site
    if(input.algorithm == ONESITE) {
      cout << "\t\toptimizing wavefunction: 1-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(lHopr[k],getfile("left-H",input.prefix,i,iroot,k));
        load(lSopr[k],getfile("left-S",input.prefix,i,iroot,k));
      }
//    eswp = optimize_onesite_merged(0,E,      OpSys,lHopr,rHopr,lSopr,rSopr,WfEnv,WfSys,0.1*T,M);
      eswp = optimize_onesite(0,E,      OpSys,lHopr,rHopr,lSopr,rSopr,WfMps,WfEnv,WfSys,0.1*T,M);
    }
    else {
      cout << "\t\toptimizing wavefunction: 2-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(lHopr[k],getfile("left-H",input.prefix,i-1,iroot,k));
        load(lSopr[k],getfile("left-S",input.prefix,i-1,iroot,k));
      }
//    eswp = optimize_twosite_merged(0,E,OpEnv,OpSys,lHopr,rHopr,lSopr,rSopr,WfEnv,WfSys,0.1*T,M);
      eswp = optimize_twosite(0,E,OpEnv,OpSys,lHopr,rHopr,lSopr,rSopr,WfMps,WfEnv,WfSys,0.1*T,M);
    }
    if(eswp < emin) emin = eswp;
    // print result
    cout << "\t\t--------------------------------------------------------------------------------" << endl;
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    cout << "\t\t--------------------------------------------------------------------------------" << endl;

    cout << "\t\tsaving operators and wavefunction for this site (sys)" << endl;
    save(WfSys[iroot],getfile("rmps",input.prefix,i,iroot));
    for(size_t k = 0; k <= iroot; ++k) {
      save(rHopr[k],getfile("right-H",input.prefix,i-1,iroot,k));
      save(rSopr[k],getfile("right-S",input.prefix,i-1,iroot,k));
    }
    // load WfSys for the next loop
    for(size_t k = 0; k < iroot; ++k) {
      load(WfSys[k],getfile("wave",input.prefix,i-1,k));
    }
    OpSys = OpEnv;
    WfSys[iroot] = WfEnv[iroot];
    save(WfSys[iroot],getfile("wave",input.prefix,i-1,iroot));
  }

  cout << "\t====================================================================================================" << endl;

  return emin;
}
