#ifndef __MPSXX_CANONICALIZE_HPP
#define __MPSXX_CANONICALIZE_HPP

#include <random>
#include <functional>

#include <btas/QSPARSE/QSTArray.h>

namespace mpsxx {

template<size_t N, class Q>
void make_noise (btas::QSTArray<double,N,Q>& wfnc, const double& noise)
{
  if(noise < 1.0e-8) return;

  std::mt19937 rgen;
  std::uniform_real_distribution<double> dist(-1.0,1.0);

  auto ds = wfnc.check_net_dshape();
  auto qs = wfnc.qshape();
  btas::QSTArray<double,N,Q> wfpt(wfnc.q(),qs,ds,std::bind(dist,rgen));

  btas::Axpy(noise,wfpt,wfnc);
}

//! Canonicalize two-site MPS
template<class Q>
void canonicalize (
        bool forward,
        btas::QSTArray<double,4,Q>& wfnc,
        btas::QSTArray<double,3,Q>& lmps,
        btas::QSTArray<double,3,Q>& rmps,
  const int& M = 0, const double& noise = 0.0)
{
  make_noise(wfnc,noise);

  btas::STArray<double,1> s;
  if(forward) {
    btas::Gesvd<double,4,3,Q,btas::LeftArrow>(wfnc,s,lmps,rmps,M);
    btas::Dimm(s,rmps);
  }
  else {
    btas::Gesvd<double,4,3,Q,btas::RightArrow>(wfnc,s,lmps,rmps,M);
    btas::Dimm(lmps,s);
  }
}

//! Canonicalize one-site MPS
template<class Q>
void canonicalize (
        bool forward,
        btas::QSTArray<double,3,Q>& wfnc,
        btas::QSTArray<double,3,Q>& lmps,
        btas::QSTArray<double,3,Q>& rmps,
  const int& M = 0, const double& noise = 0.0)
{
  if(noise >= 1.0e-8) {
    btas::QSTArray<double,4,Q> temp;
    if(forward) {
      btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,wfnc,rmps,1.0,temp);
    }
    else {
      btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lmps,wfnc,1.0,temp);
    }
    return canonicalize(forward,temp,lmps,rmps,M,noise);
  }

  btas::STArray<double,1> s;
  btas::QSTArray<double,2,Q> g;
  if(forward) {
    btas::Gesvd<double,3,3,Q,btas::LeftArrow>(wfnc,s,lmps,g,M);
    btas::Dimm(s,g);
    btas::QSTArray<double,3,Q> temp;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,g,rmps,1.0,temp);
    rmps = temp;
  }
  else {
    btas::Gesvd<double,3,2,Q,btas::RightArrow>(wfnc,s,g,rmps,M);
    btas::Dimm(g,s);
    btas::QSTArray<double,3,Q> temp;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lmps,g,1.0,temp);
    lmps = temp;
  }
}

//! Canonicalize merged MPS
template<class Q>
void canonicalize (
        bool forward,
        btas::QSTArray<double,2,Q>& wfnc,
        btas::QSTArray<double,2,Q>& lmat,
        btas::QSTArray<double,2,Q>& rmat,
  const int& M = 0, const double& noise = 0.0)
{
  make_noise(wfnc,noise);

  btas::STArray<double,1> s;
  if(forward) {
    btas::Gesvd<double,2,2,Q,btas::LeftArrow>(wfnc,s,lmat,rmat,M);
    btas::Dimm(s,rmat);
  }
  else {
    btas::Gesvd<double,2,2,Q,btas::RightArrow>(wfnc,s,lmat,rmat,M);
    btas::Dimm(lmat,s);
  }
}

} // namespace mpsxx

#endif // __MPSXX_CANONICALIZE_HPP
