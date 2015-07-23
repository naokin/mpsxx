#ifndef __MPSXX_DMRG_DENSE_INITIALIZE_HPP
#define __MPSXX_DMRG_DENSE_INITIALIZE_HPP

#include <random>

#include <boost/bind.hpp>

#include "MPX.h"
#include "driver.hpp"

namespace mpsxx {

template<typename T>
void gen_random_mpss (size_t N, MPSs<T>& mpss, size_t M = 1)
{
  mpss.resize(N);

  std::uniform_real_distribution<double> dist(-1.0,1.0);
  std::mt19937 rgen;

  size_t Mx = 1;

  for(size_t i = 0; i < N-1; ++i) {
    mpss[i].resize(Mx,4,M);
    mpss[i].generate(boost::bind(dist,rgen));
    Mx = M;
  }

  mpss[N-1].resize(Mx,4,1);
  mpss[N-1].generate(boost::bind(dist,rgen));
}

/// canonicalize & renormalize
template<typename T>
void initialize (const MPOs<T>& mpos, MPSs<T>& mpss, std::vector<BLOCK<T>>& lops, std::vector<BLOCK<T>>& rops, size_t M = 1)
{
  size_t N = mpos.size();

  gen_random_mpss(N,mpss,M);

  rops.resize(N);
  rops[N-1].resize(1,1,1);
  rops[N-1].fill(1.0);

  for(size_t i = N-1; i > 0; --i) {
    btas::Normalize(mpss[i]);
    btas::TArray<T,3> wfnc(mpss[i]);
    btas::TArray<T,2> rest;
    canonicalize(0,wfnc,mpss[i],rest,M);
    wfnc.clear();
    btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,mpss[i-1],rest,1.0,wfnc);
    mpss[i-1] = wfnc;
    renormalize (0,mpos[i],rops[i],mpss[i],mpss[i],rops[i-1]);
  }

  btas::Normalize(mpss[0]);

  lops.resize(N);
  lops[0].resize(1, 1, 1);
  lops[0].fill(1.0);
}

} // namespace mpsxx

#endif // __MPSXX_DMRG_DENSE_INITIALIZE_HPP
