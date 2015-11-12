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
  using std::endl;
  using std::flush;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  Communicator world;

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

  pout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  pout << "\t\t\tFORWARD SWEEP" << endl;
  pout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  load(OpSys,getfile("mpo",input.prefix,0));
  for(size_t k = 0; k <= iroot; ++k) {
    if(world.rank() == 0) {
      load(WfSys[k],getfile("wave",input.prefix,0,k));
      load(lSopr[k],getfile("left-S",input.prefix,0,iroot,k));
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,WfSys[k],0);
    boost::mpi::broadcast(world,lSopr[k],0);
#endif
    load(lHopr[k],getfile("left-H",input.prefix,0,iroot,k));
  }

  for(size_t i = 0; i < N-1; ++i) {

    pout << "\t====================================================================================================" << endl;
    pout << "\t\tSITE [ " << setw(3) << i << " ] " << endl;
    pout << "\t----------------------------------------------------------------------------------------------------" << endl;

    pout << "\t\tloading operators and wavefunction of next site (env)" << endl;
    load(OpEnv,getfile("mpo",input.prefix,i+1));
    for(size_t k = 0; k <= iroot; ++k) {
      if(world.rank() == 0) {
        load(WfEnv[k],getfile("rmps",input.prefix,i+1,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,WfEnv[k],0);
#endif
    }
    for(size_t k = 0; k < iroot; ++k) {
      if(world.rank() == 0) {
        load(WfMps[k],getfile("lmps",input.prefix,i,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,WfMps[k],0);
#endif
    }

    double eswp;

    // Optimize for each site
    if(input.algorithm == ONESITE) {
      pout << "\t\toptimizing wavefunction: 1-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(rHopr[k],getfile("right-H",input.prefix,i,iroot,k));
        if(world.rank() == 0) {
          load(rSopr[k],getfile("right-S",input.prefix,i,iroot,k));
        }
#ifndef _SERIAL
        boost::mpi::broadcast(world,rSopr[k],0);
#endif
      }
//    eswp = optimize_onesite_merged(1,E,OpSys,      lHopr,rHopr,lSopr,rSopr,WfSys,WfEnv,0.1*T,M,input.noise);
      eswp = optimize_onesite(1,E,OpSys,      lHopr,rHopr,lSopr,rSopr,WfMps,WfSys,WfEnv,0.1*T,M,input.noise);
    }
    else {
      pout << "\t\toptimizing wavefunction: 2-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(rHopr[k],getfile("right-H",input.prefix,i+1,iroot,k));
        if(world.rank() == 0) {
          load(rSopr[k],getfile("right-S",input.prefix,i+1,iroot,k));
        }
#ifndef _SERIAL
        boost::mpi::broadcast(world,rSopr[k],0);
#endif
      }
//    eswp = optimize_twosite_merged(1,E,OpSys,OpEnv,lHopr,rHopr,lSopr,rSopr,WfSys,WfEnv,0.1*T,M,input.noise);
      eswp = optimize_twosite(1,E,OpSys,OpEnv,lHopr,rHopr,lSopr,rSopr,WfMps,WfSys,WfEnv,0.1*T,M,input.noise);
    }
    if(eswp < emin) emin = eswp;
    // print result
    pout << "\t\t--------------------------------------------------------------------------------" << endl;
    pout.precision(16);
    pout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    pout << "\t\t--------------------------------------------------------------------------------" << endl;

    pout << "\t\tsaving operators and wavefunction of this site (sys)" << endl;
    if(world.rank() == 0) {
      save(WfSys[iroot],getfile("lmps",input.prefix,i,iroot));
    }
    for(size_t k = 0; k <= iroot; ++k) {
      save(lHopr[k],getfile("left-H",input.prefix,i+1,iroot,k));
      if(world.rank() == 0) {
        save(lSopr[k],getfile("left-S",input.prefix,i+1,iroot,k));
      }
    }
    // load WfSys for the next loop
    for(size_t k = 0; k < iroot; ++k) {
      if(world.rank() == 0) {
        load(WfSys[k],getfile("wave",input.prefix,i+1,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,WfSys[k],0);
#endif
    }
    OpSys = OpEnv;
    WfSys[iroot] = WfEnv[iroot];
    if(world.rank() == 0) {
      save(WfSys[iroot],getfile("wave",input.prefix,i+1,iroot));
    }
  }

  pout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  pout << "\t\t\tBACKWARD SWEEP" << endl;
  pout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  for(size_t k = 0; k <= iroot; ++k) {
    load(rHopr[k],getfile("right-H",input.prefix,N-1,iroot,k));
    if(world.rank() == 0) {
      load(rSopr[k],getfile("right-S",input.prefix,N-1,iroot,k));
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,rSopr[k],0);
#endif
  }

  for(size_t i = N-1; i > 0; --i) {

    pout << "\t====================================================================================================" << endl;
    pout << "\t\tSITE [ " << setw(3) << i << " ] " << endl;
    pout << "\t----------------------------------------------------------------------------------------------------" << endl;

    pout << "\t\tloading operators and wavefunction of next site (env)" << endl;
    load(OpEnv,getfile("mpo",input.prefix,i-1));
    for(size_t k = 0; k <= iroot; ++k) {
      if(world.rank() == 0) {
        load(WfEnv[k],getfile("lmps",input.prefix,i-1,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,WfEnv[k],0);
#endif
    }
    for(size_t k = 0; k < iroot; ++k) {
      if(world.rank() == 0) {
        load(WfMps[k],getfile("rmps",input.prefix,i,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,WfMps[k],0);
#endif
    }

    double eswp;

    // Optimize for each site
    if(input.algorithm == ONESITE) {
      pout << "\t\toptimizing wavefunction: 1-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(lHopr[k],getfile("left-H",input.prefix,i,iroot,k));
        if(world.rank() == 0) {
          load(lSopr[k],getfile("left-S",input.prefix,i,iroot,k));
        }
#ifndef _SERIAL
        boost::mpi::broadcast(world,rSopr[k],0);
#endif
      }
//    eswp = optimize_onesite_merged(0,E,      OpSys,lHopr,rHopr,lSopr,rSopr,WfEnv,WfSys,0.1*T,M,input.noise);
      eswp = optimize_onesite(0,E,      OpSys,lHopr,rHopr,lSopr,rSopr,WfMps,WfEnv,WfSys,0.1*T,M,input.noise);
    }
    else {
      pout << "\t\toptimizing wavefunction: 2-site algorithm " << endl;
      for(size_t k = 0; k <= iroot; ++k) {
        load(lHopr[k],getfile("left-H",input.prefix,i-1,iroot,k));
        if(world.rank() == 0) {
          load(lSopr[k],getfile("left-S",input.prefix,i-1,iroot,k));
        }
#ifndef _SERIAL
        boost::mpi::broadcast(world,rSopr[k],0);
#endif
      }
//    eswp = optimize_twosite_merged(0,E,OpEnv,OpSys,lHopr,rHopr,lSopr,rSopr,WfEnv,WfSys,0.1*T,M,input.noise);
      eswp = optimize_twosite(0,E,OpEnv,OpSys,lHopr,rHopr,lSopr,rSopr,WfMps,WfEnv,WfSys,0.1*T,M,input.noise);
    }
    if(eswp < emin) emin = eswp;
    // print result
    pout << "\t\t--------------------------------------------------------------------------------" << endl;
    pout.precision(16);
    pout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    pout << "\t\t--------------------------------------------------------------------------------" << endl;

    pout << "\t\tsaving operators and wavefunction for this site (sys)" << endl;
    if(world.rank() == 0) {
      save(WfSys[iroot],getfile("rmps",input.prefix,i,iroot));
    }
    for(size_t k = 0; k <= iroot; ++k) {
      save(rHopr[k],getfile("right-H",input.prefix,i-1,iroot,k));
      if(world.rank() == 0) {
        save(rSopr[k],getfile("right-S",input.prefix,i-1,iroot,k));
      }
    }
    // load WfSys for the next loop
    for(size_t k = 0; k < iroot; ++k) {
      if(world.rank() == 0) {
        load(WfSys[k],getfile("wave",input.prefix,i-1,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,WfSys[k],0);
#endif
    }
    OpSys = OpEnv;
    WfSys[iroot] = WfEnv[iroot];
    if(world.rank() == 0) {
      save(WfSys[iroot],getfile("wave",input.prefix,i-1,iroot));
    }
  }

  pout << "\t====================================================================================================" << endl;

  return emin;
}
