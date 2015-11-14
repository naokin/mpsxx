#include <iostream>
#include <iomanip>
#include <random>

#include "make_random_mpss.h"

#include "mpidefs.h"

#include "mpsite.h"
#include "canonicalize.hpp"
#include "renormalize.hpp"
#include "fileio.h"
#include "fermion.h"

/// Generate quantum numbers for each boundary from MPO
///
/// \param sites quantum numbers for each site (physical index)
/// \param qt total quantum number of lattice
/// \param max_quantum_blocks maximum number of quantum blocks on each boundary (0: no limitation by default)
std::vector<btas::Qshapes<fermion>>
mpsxx::generate_quantum_states (
  const std::vector<btas::Qshapes<fermion>>& sites,
  const fermion& qt, size_t max_quantum_blocks)
{
  size_t N = sites.size();

  // zero quantum number
  btas::Qshapes<fermion> qz(1, fermion::zero());

  // boundary quantum numbers
  std::vector<btas::Qshapes<fermion>> qb(N);

  // generate quantum number blocks for MPS
  qb[0] = sites[0];
  for(size_t i = 1; i < N-1; ++i) {
    btas::Qshapes<fermion> qb_bare = qb[i-1] & sites[i]; // get unique elements of { q(l) x q(n) }
    qb[i].clear();
    qb[i].reserve(qb_bare.size());
    for(size_t j = 0; j < qb_bare.size(); ++j)
      if(qb_bare[j].p() <= qt.p()) qb[i].push_back(qb_bare[j]);
  }
  qb[N-1] = btas::Qshapes<fermion>(1, qt); // used for checking quantum states

  // reduce zero quantum blocks
  for(size_t i = N-1; i > 0; --i) {
    const btas::Qshapes<fermion>& qn = sites[i];
          btas::Qshapes<fermion>& ql = qb[i-1];
          btas::Qshapes<fermion>& qr = qb[i];

    // check non-zero for each ql index
    auto lt = ql.begin();
    while(lt != ql.end()) {
      bool non_zero = false;
      for(size_t p = 0; p < qn.size(); ++p) {
        for(size_t r = 0; r < qr.size(); ++r) {
          non_zero |= (qr[r] == (qn[p] * (*lt)));
        }
      }
      if(non_zero)
        ++lt;
      else
        ql.erase(lt);
    }
    assert(ql.size() > 0);
  }
  for(size_t i = 0; i < N-1; ++i) {
    const btas::Qshapes<fermion>& qn = sites[i];
          btas::Qshapes<fermion>& ql = qb[i];
          btas::Qshapes<fermion>& qr = qb[i+1];

    // further reduction
    if(max_quantum_blocks > 0 && ql.size() > max_quantum_blocks) {
      size_t offs = (ql.size()-max_quantum_blocks)/2;
      ql = btas::Qshapes<fermion>(ql.begin()+offs, ql.begin()+offs+max_quantum_blocks);

      // check non-zero for each qr index
      auto rt = qr.begin();
      while(rt != qr.end()) {
        bool non_zero = false;
        for(size_t l = 0; l < ql.size(); ++l) {
          for(size_t p = 0; p < qn.size(); ++p) {
            non_zero |= (*rt == (ql[l] * qn[p]));
          }
        }
        if(non_zero)
          ++rt;
        else
          qr.erase(rt);
      }
      assert(qr.size() > 0);
    }
  }
  return qb;
}

/// Generate random MPS for a specified root, \c iroot
void mpsxx::make_random_mpss (const mpsxx::DMRGInput& input, const size_t& iroot)
{
  using std::endl;
  using std::setw;

  Communicator world;

  size_t N = input.N_sites;
  size_t K = input.N_roots;
  int    M = input.N_max_states;

  std::mt19937 rgen;
  std::uniform_real_distribution<double> dist(-1.0,1.0);

  fermion qt(input.N_elecs,input.N_spins);

  pout << "\t====================================================================================================" << endl;
  pout << "\t\t\tGenerate & initialize MPS: number of sites = " << setw(3) << N << endl;
  pout << "\t====================================================================================================" << endl;

  btas::QSTArray<double,3,fermion> wfnc;

  btas::Qshapes<fermion> qz(1,fermion::zero());

  if(world.rank() == 0) {

  btas::Qshapes<fermion> ql(1,fermion::zero());
  btas::Dshapes          dl(ql.size(),1);

  btas::Qshapes<fermion> qr;
  btas::Dshapes          dr;

  std::vector<btas::Qshapes<fermion>> qn(N);
  std::vector<btas::Dshapes>          dn(N);

  for(size_t i = 0; i < N; ++i) {
    qn[i] = MpSite<fermion>::quanta();
    dn[i] = btas::Dshapes(qn[i].size(),1);
  }

  pout << "\t\t\tGenerating quantum states for each boundary " << endl;
  std::vector<btas::Qshapes<fermion>> qb = generate_quantum_states(qn,qt,M);
  pout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  btas::TVector<btas::Qshapes<fermion>,3> q_shape;
  btas::TVector<btas::Dshapes,         3> d_shape;

  for(size_t i = 0; i < N-1; ++i) {
    pout << "\t\t\tGenerating   MPS for site [ " << setw(3) << i << " ] " << endl;
    pout << "\t----------------------------------------------------------------------------------------------------" << endl;
    qr = qb[i];
    dr = btas::Dshapes(qr.size(),1);

    q_shape = btas::make_array(ql,qn[i],-qr);
    d_shape = btas::make_array(dl,dn[i], dr);

    wfnc.resize(fermion::zero(),q_shape,d_shape,dist(rgen));
    btas::Normalize(wfnc);
    save(wfnc,getfile("wave",input.prefix,i,iroot));
    wfnc.clear();

    ql = qr;
    dl = dr;
  }

  pout << "\t\t\tGenerating   MPS for site [ " << setw(3) << N-1 << " ] " << endl;
  pout << "\t====================================================================================================" << endl;

  qr = qz; // qb[N-1] is not used
  dr = btas::Dshapes(qr.size(),1);

  q_shape = btas::make_array(ql,qn[N-1],-qr);
  d_shape = btas::make_array(dl,dn[N-1], dr);

  wfnc.resize(qt,q_shape,d_shape,dist(rgen)); // set qt as a total quantum number of array
  btas::Normalize(wfnc);
  save(wfnc,getfile("wave",input.prefix,N-1,iroot));

  } // if(world.rank() == 0)

  // peak group size
  size_t ngroup = 0;
  {
    std::vector<btas::QSTArray<double,4,fermion>> dummy;
    load(dummy,getfile("mpo",input.prefix,0));
    ngroup = dummy.size();
  }

  // Initial canonicalization
  std::vector<std::vector<btas::QSTArray<double,3,fermion>>> rH0(iroot+1);
  std::vector<btas::QSTArray<double,2,fermion>> rS0(iroot+1);
  for(size_t k = 0; k <= iroot; ++k) {
    rH0[k].resize(ngroup);
    for(size_t g = 0; g < ngroup; ++g) {
      rH0[k][g].resize(fermion::zero(),btas::make_array(qz,qz,qz));
      rH0[k][g].insert(btas::shape(0,0,0),btas::TArray<double,3>(1,1,1));
      rH0[k][g].fill(1.0);
    }
    save(rH0[k],getfile("right-H",input.prefix,N-1,iroot,k));
    if(world.rank() == 0) {
      rS0[k].resize(fermion::zero(),btas::make_array(qz,qz));
      rS0[k].insert(btas::shape(0,0),btas::TArray<double,2>(1,1));
      rS0[k].fill(1.0);
      save(rS0[k],getfile("right-S",input.prefix,N-1,iroot,k));
    }
  }

  for(size_t i = N-1; i > 0; --i) {
    pout << "\t\t\tInitializing MPS for site [ " << setw(3) << i << " ] " << endl;
    pout << "\t----------------------------------------------------------------------------------------------------" << endl;
    btas::QSTArray<double,3,fermion> rmps;
    btas::QSTArray<double,3,fermion> lmps;
    if(world.rank() == 0) {
      load(lmps,getfile("wave",input.prefix,i-1,iroot));

      canonicalize(0,wfnc,lmps,rmps,M);
      wfnc = lmps;
      btas::Normalize(wfnc);

      save(wfnc,getfile("wave",input.prefix,i-1,iroot));
      save(rmps,getfile("rmps",input.prefix,i,iroot));
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,wfnc,0);
    boost::mpi::broadcast(world,rmps,0);
#endif

    std::vector<btas::QSTArray<double,4,fermion>> mpo;
    load(mpo,getfile("mpo",input.prefix,i));
    assert(mpo.size() == ngroup);

    for(size_t k = 0; k <= iroot; ++k) {
      btas::QSTArray<double,3,fermion> kmps;
      if(world.rank() == 0) {
        load(kmps,getfile("rmps",input.prefix,i,k));
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,kmps,0);
#endif

      std::vector<btas::QSTArray<double,3,fermion>> rH1(ngroup);
      for(size_t g = 0; g < ngroup; ++g) {
        renormalize(0,mpo[g],rH0[k][g],rmps,kmps,rH1[g]);
      }
      save(rH1,getfile("right-H",input.prefix,i-1,iroot,k));
      rH0[k] = rH1;

      if(world.rank() == 0) {
        btas::QSTArray<double,2,fermion> rS1;
        renormalize(0,rS0[k],rmps,kmps,rS1);
        save(rS1,getfile("right-S",input.prefix,i-1,iroot,k));
        rS0[k] = rS1;
      }
    }
  }
  pout << "\t\t\tInitializing MPS for site [ " << setw(3) << 0 << " ] " << endl;
  pout << "\t====================================================================================================" << endl;

  std::vector<std::vector<btas::QSTArray<double,3,fermion>>> lH0(iroot+1);
  std::vector<btas::QSTArray<double,2,fermion>> lS0(iroot+1);
  for(size_t k = 0; k <= iroot; ++k) {
    lH0[k].resize(ngroup);
    for(size_t g = 0; g < ngroup; ++g) {
      lH0[k][g].resize(fermion::zero(),btas::make_array(qz,qz,qz));
      lH0[k][g].insert(btas::shape(0,0,0),btas::TArray<double,3>(1,1,1));
      lH0[k][g].fill(1.0);
    }
    save(lH0[k],getfile("left-H",input.prefix,0,iroot,k));
    if(world.rank() == 0) {
      lS0[k].resize(fermion::zero(),btas::make_array(qz,qz));
      lS0[k].insert(btas::shape(0,0),btas::TArray<double,2>(1,1));
      lS0[k].fill(1.0);
      save(lS0[k],getfile("left-S",input.prefix,0,iroot,k));
    }
  }
}
