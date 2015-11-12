#ifndef __MPSXX_EXPECTATION_HPP
#define __MPSXX_EXPECTATION_HPP

#include <iostream>
#include <iomanip>

#include <legacy/QSPARSE/QSTArray.h>

#include "input.h"
#include "renormalize.hpp"
#include "fileio.h"

namespace mpsxx {

/// Compute overlap b/w two different states <i|j>
template<class Q>
double overlap (const DMRGInput& input, const size_t& iroot = 0, const size_t& jroot = 0)
{
  using std::cout;
  using std::endl;
  using std::setw;

  size_t N = input.N_sites;
  size_t K = input.N_roots;

  assert(iroot < K);
  assert(jroot < K);

  btas::QSTArray<double,3,Q> mpsi;
  btas::QSTArray<double,3,Q> mpsj;
  btas::QSTArray<double,2,Q> lSij;

  btas::Qshapes<Q> qz(1,Q::zero());
  lSij.resize(Q::zero(),btas::make_array(qz,qz));
  lSij.insert(btas::shape(0,0),btas::TArray<double,2>(1,1));
  lSij.fill(1.0);

  for(size_t i = 0; i < N-1; ++i) {
    load(mpsi,getfile("lmps",input.prefix,i,iroot));
    load(mpsj,getfile("lmps",input.prefix,i,jroot));
    btas::QSTArray<double,2,Q> lSijTmp;
    renormalize(1,lSij,mpsi,mpsj,lSijTmp);
    lSij = lSijTmp;
  }
    load(mpsi,getfile("wave",input.prefix,N-1,iroot));
    load(mpsj,getfile("wave",input.prefix,N-1,jroot));
    btas::QSTArray<double,2,Q> tSij;
    renormalize(1,lSij,mpsi,mpsj,tSij);

  return tSij.find(0)->second->data()[0];
}

/// Compute expectation value <i|O|j>
template<class Q>
double expectation (const DMRGInput& input, const size_t& iroot = 0, const size_t& jroot = 0)
{
  using std::cout;
  using std::endl;
  using std::setw;

  size_t N = input.N_sites;
  size_t K = input.N_roots;

  assert(iroot < K);
  assert(jroot < K);

  btas::QSTArray<double,4,Q> mpo;
  btas::QSTArray<double,3,Q> mpsi;
  btas::QSTArray<double,3,Q> mpsj;
  btas::QSTArray<double,3,Q> lHij;

  btas::Qshapes<Q> qz(1,Q::zero());
  lHij.resize(Q::zero(),btas::make_array(qz,qz,qz));
  lHij.insert(btas::shape(0,0,0),btas::TArray<double,3>(1,1,1));
  lHij.fill(1.0);

  for(size_t i = 0; i < N-1; ++i) {
    load(mpo, getfile("mpo", input.prefix,i));
    load(mpsi,getfile("lmps",input.prefix,i,iroot));
    load(mpsj,getfile("lmps",input.prefix,i,jroot));
    btas::QSTArray<double,3,Q> lHijTmp;
    renormalize(1,mpo,lHij,mpsi,mpsj,lHijTmp);
    lHij = lHijTmp;
  }
    load(mpo, getfile("mpo", input.prefix,N-1));
    load(mpsi,getfile("wave",input.prefix,N-1,iroot));
    load(mpsj,getfile("wave",input.prefix,N-1,jroot));
    btas::QSTArray<double,3,Q> tHij;
    renormalize(1,mpo,lHij,mpsi,mpsj,tHij);

  return tHij.find(0)->second->data()[0];

//btas::QSTArray<double,3,Q> rHij;

//rHij.resize(Q::zero(),btas::make_array(qz,qz,qz));
//rHij.insert(btas::shape(0,0,0),btas::TArray<double,3>(1,1,1));
//rHij.fill(1.0);

//for(size_t i = N-1; i > 0; --i) {
//  load(mpo, getfile("mpo", input.prefix,i));
//  load(mpsi,getfile("rmps",input.prefix,i,iroot));
//  load(mpsj,getfile("rmps",input.prefix,i,jroot));
//  btas::QSTArray<double,3,Q> rHijTmp;
//  renormalize(0,mpo,rHij,mpsi,mpsj,rHijTmp);
//  rHij = rHijTmp;
//}
//  load(mpo, getfile("mpo", input.prefix,0));
//  load(mpsi,getfile("wave",input.prefix,0,iroot));
//  load(mpsj,getfile("wave",input.prefix,0,jroot));
//  btas::QSTArray<double,3,Q> uHij;
//  renormalize(0,mpo,rHij,mpsi,mpsj,uHij);

//return uHij.find(0)->second->data()[0];
}

} // namespace mpsxx

#endif // __MPSXX_EXPECTATION_HPP
