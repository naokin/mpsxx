#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <boost/filesystem.hpp>

#include "make_qc_mpos.h"

#include "mpidefs.h"

#include "fileio.h"
#include "gen_qc_operators.h"
#include "compress_qc_mpos.h"

int mpsxx::get_group (const size_t& Ndiv, const int& Norbs, const int& i, const int& j)
{
  int ix = i;
  int jx = j;
  if(ix > jx) std::swap(ix,jx);

  return (jx*(jx+1)/2+ix)%Ndiv;
}

int mpsxx::get_group (const size_t& Ndiv, const int& Norbs, const int& i, const int& j, const int& k, const int& l)
{
  std::vector<int> index = { i, j, k, l };
  std::sort(index.begin(),index.end());
  const int& ix = index[0];
  const int& jx = index[1];
  const int& kx = Norbs-index[2];
  const int& lx = Norbs-index[3];

  int ij = jx*(jx+1)/2+ix;
  int kl = kx*(kx+1)/2+lx;

  if(ij < kl)
    return ij%Ndiv;
  else
    return kl%Ndiv;
}

void mpsxx::make_qc_mpos (
  const size_t& Ndivs,
  const size_t& Norbs,
  const double& Ecore,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint,
  const std::string& opname,
  const std::string& prefix,
  const bool& enable_swap_sweep,
  const bool& do_compress)
{
  std::vector<int> group;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  btas::TArray<int,2> grp1e(Norbs,Norbs);
  grp1e.fill(-1);
  btas::TArray<int,4> grp2e(Norbs,Norbs,Norbs,Norbs);
  grp2e.fill(-1);

  if(world.rank() == 0) {

    std::set<int> gset;

    for(int i = 0; i < Norbs; ++i)
      for(int j = 0; j < Norbs; ++j) {
        if(fabs(oneint(i,j)) >= 1.0e-16) {
          int g = get_group(Ndivs,Norbs,i,j);
          if(g >= 0) gset.insert(g);
          grp1e(i,j) = g;
        }
      }
    for(int i = 0; i < Norbs; ++i)
      for(int k = 0; k < Norbs; ++k)
        for(int j = 0; j < Norbs; ++j)
          for(int l = 0; l < Norbs; ++l) {
            if(fabs(twoint(i,k,j,l)) >= 1.0e-16) {
              int g = get_group(Ndivs,Norbs,i,k,j,l);
              if(g >= 0) gset.insert(g);
              grp2e(i,k,j,l) = g;
            }
          }
    group.assign(gset.begin(),gset.end());
  }
#ifndef _SERIAL
  boost::mpi::broadcast(world,group,0);
  boost::mpi::broadcast(world,grp1e,0);
  boost::mpi::broadcast(world,grp2e,0);
#endif

  size_t g_local = 0;

  for(size_t g = 0; g < group.size(); ++g) {

    double e = 0.0; if(g == 0) e = Ecore;

    if(g%nproc == iproc) {
      btas::TArray<double,2> tmp1e(Norbs,Norbs);
      tmp1e.fill(0.0);
      btas::TArray<double,4> tmp2e(Norbs,Norbs,Norbs,Norbs);
      tmp2e.fill(0.0);
      for(int i = 0; i < Norbs; ++i)
        for(int j = 0; j < Norbs; ++j)
          if(grp1e(i,j) == group[g])
            tmp1e(i,j) = oneint(i,j);
      for(int i = 0; i < Norbs; ++i)
        for(int k = 0; k < Norbs; ++k)
          for(int j = 0; j < Norbs; ++j)
            for(int l = 0; l < Norbs; ++l)
              if(grp2e(i,k,j,l) == group[g])
                tmp2e(i,k,j,l) = twoint(i,k,j,l);

      std::ostringstream oss;
      oss << "mposcr-" << g_local;
      gen_qc_operators(Norbs,e,tmp1e,tmp2e,oss.str(),prefix,enable_swap_sweep);

      if(do_compress) compress_qc_mpos(Norbs,oss.str(),prefix);

      ++g_local;
    }
  }

  for(size_t s = 0; s < Norbs; ++s) {
    std::vector<btas::QSTArray<double,4,fermion>> mpo(g_local);
    for(size_t g = 0; g < g_local; ++g) {
      std::ostringstream oss;
      oss << "mposcr-" << g;
      load(mpo[g],getfile(oss.str(),prefix,s));

      // remove mposcr....
      boost::filesystem::path path_to_scr(getfile(oss.str(),prefix,s));
      assert(boost::filesystem::exists(path_to_scr));
      boost::filesystem::remove(path_to_scr);
    }
    save(mpo,getfile(opname,prefix,s));
  }
}
