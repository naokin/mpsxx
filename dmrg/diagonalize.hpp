#ifndef __MPSXX_DAVIDSON_DIAGONALIZE_HPP
#define __MPSXX_DAVIDSON_DIAGONALIZE_HPP

#include <iostream>
#include <iomanip>

#include <vector>
#include <cmath>
#include <functional>

#include <legacy/QSPARSE/QSTArray.h>

#include "mpidefs.h"
#include "compute_sigma_vector.hpp"

namespace mpsxx {
namespace Davidson {

template<size_t N, class Q> struct Functor;

template<class Q>
struct Functor<3,Q> {

  typedef btas::QSTArray<double,3,Q> MPS_type;

  typedef btas::QSTArray<double,4,Q> MPO_type;

  typedef btas::QSTArray<double,3,Q> HOp_type;

  typedef btas::QSTArray<double,2,Q> SOp_type;

  Functor (
    const std::vector<double>& E,
    const MPO_type& mpo,
    const std::vector<HOp_type>& lHopr,
    const std::vector<HOp_type>& rHopr,
    const std::vector<SOp_type>& lSopr,
    const std::vector<SOp_type>& rSopr,
    const std::vector<MPS_type>& wfn)
  {
    N_ = wfn.size();

    if(N_ > 1) {
      E_.resize(N_-1);
      H_.resize(N_-1);
      S_.resize(N_-1);
    }

    pMPO_ = &mpo;
    pLOp_ = &lHopr[N_-1];
    pROp_ = &rHopr[N_-1];

    for(size_t k = 0; k < N_-1; ++k) {
      compute_sigma_vector(mpo,lHopr[k],rHopr[k],wfn[k],H_[k]);
      compute_sigma_vector(    lSopr[k],rSopr[k],wfn[k],S_[k]);
      E_[k] = E[k];
    }
  }

  // Hk = (1-Pi)H(1-Pi) = H - PiH - HPi + PiHPi
  // thus,
  // Hk|k> = H|k> - |i><i|H|k> - H|i><i|k> + |i><i|H|i><i|k>
  double operator() (const MPS_type& c, MPS_type& v) const
  {
    Communicator world;
    compute_sigma_vector(*pMPO_,*pLOp_,*pROp_,c,v);
    for(size_t k = 0; k < N_-1; ++k) {
      double Hk = btas::Dotc(H_[k],c);
      double Sk = btas::Dotc(S_[k],c);
      btas::Axpy(-Hk,S_[k],v);
      btas::Axpy(-Sk,H_[k],v);
      if(world.rank() == 0) {
        btas::Axpy(E_[k]*Sk,S_[k],v);
      }
    }
  }

  void orthonormalize (MPS_type& c) const
  {
    for(size_t k = 0; k < N_-1; ++k) {
      double Sk = btas::Dotc(S_[k],c);
      btas::Axpy(-Sk,S_[k],c);
    }
    btas::Normalize(c);
  }

  const MPO_type* pMPO_;

  const HOp_type* pLOp_;

  const HOp_type* pROp_;

  size_t N_;

  std::vector<double> E_;

  std::vector<MPS_type> H_;

  std::vector<MPS_type> S_;

};

template<class Q>
struct Functor<4,Q> {

  typedef btas::QSTArray<double,4,Q> MPS_type;

  typedef btas::QSTArray<double,4,Q> MPO_type;

  typedef btas::QSTArray<double,3,Q> HOp_type;

  typedef btas::QSTArray<double,2,Q> SOp_type;

  Functor (
    const std::vector<double>& E,
    const MPO_type& lmpo,
    const MPO_type& rmpo,
    const std::vector<HOp_type>& lHopr,
    const std::vector<HOp_type>& rHopr,
    const std::vector<SOp_type>& lSopr,
    const std::vector<SOp_type>& rSopr,
    const std::vector<MPS_type>& wfn)
  {
    N_ = wfn.size();

    if(N_ > 1) {
      E_.resize(N_-1);
      H_.resize(N_-1);
      S_.resize(N_-1);
    }

    pLMPO_ = &lmpo;
    pRMPO_ = &rmpo;
    pLOp_ = &lHopr[N_-1];
    pROp_ = &rHopr[N_-1];

    for(size_t k = 0; k < N_-1; ++k) {
      compute_sigma_vector(lmpo,rmpo,lHopr[k],rHopr[k],wfn[k],H_[k]);
      compute_sigma_vector(          lSopr[k],rSopr[k],wfn[k],S_[k]);
      E_[k] = E[k];
    }
  }

  // Hk = (1-Pi)H(1-Pi) = H - PiH - HPi + PiHPi
  // thus,
  // Hk|k> = H|k> - |i><i|H|k> - H|i><i|k> + |i><i|H|i><i|k>
  double operator() (const MPS_type& c, MPS_type& v) const
  {
    Communicator world;
    compute_sigma_vector(*pLMPO_,*pRMPO_,*pLOp_,*pROp_,c,v);
    for(size_t k = 0; k < N_-1; ++k) {
      double Hk = btas::Dotc(H_[k],c);
      double Sk = btas::Dotc(S_[k],c);
      btas::Axpy(-Hk,S_[k],v);
      btas::Axpy(-Sk,H_[k],v);
      if(world.rank() == 0) {
        btas::Axpy(E_[k]*Sk,S_[k],v);
      }
    }
  }

  void orthonormalize (MPS_type& c) const
  {
    for(size_t k = 0; k < N_-1; ++k) {
      double Sk = btas::Dotc(S_[k],c);
      btas::Axpy(-Sk,S_[k],c);
    }
    btas::Normalize(c);
  }

  const MPO_type* pLMPO_;

  const MPO_type* pRMPO_;

  const HOp_type* pLOp_;

  const HOp_type* pROp_;

  size_t N_;

  std::vector<double> E_;

  std::vector<MPS_type> H_;

  std::vector<MPS_type> S_;

};

/// For merged MPS
template<class Q>
struct Functor<2,Q> {

  typedef btas::QSTArray<double,2,Q> MPS_type;

  typedef btas::QSTArray<double,3,Q> HOp_type;

  typedef btas::QSTArray<double,2,Q> SOp_type;

  Functor (
    const std::vector<double>& E,
    const std::vector<HOp_type>& lHopr,
    const std::vector<HOp_type>& rHopr,
    const std::vector<SOp_type>& lSopr,
    const std::vector<SOp_type>& rSopr,
    const std::vector<MPS_type>& wfn)
  {
    N_ = wfn.size();
    if(N_ > 1) {
      E_.resize(N_-1);
      H_.resize(N_-1);
      S_.resize(N_-1);
    }

    pLOp_ = &lHopr[N_-1];
    pROp_ = &rHopr[N_-1];

    for(size_t k = 0; k < N_-1; ++k) {
      compute_sigma_vector(lHopr[k],rHopr[k],wfn[k],H_[k]);
      compute_sigma_vector(lSopr[k],rSopr[k],wfn[k],S_[k]);
      E_[k] = E[k];
    }
  }

  // Hk = (1-Pi)H(1-Pi) = H - PiH - HPi + PiHPi
  // thus,
  // Hk|k> = H|k> - |i><i|H|k> - H|i><i|k> + |i><i|H|i><i|k>
  double operator() (const MPS_type& c, MPS_type& v) const
  {
    Communicator world;
    compute_sigma_vector(*pLOp_,*pROp_,c,v);
    for(size_t k = 0; k < N_-1; ++k) {
      double Hk = btas::Dotc(H_[k],c);
      double Sk = btas::Dotc(S_[k],c);
      btas::Axpy(-Hk,S_[k],v);
      btas::Axpy(-Sk,H_[k],v);
      if(world.rank() == 0) {
        btas::Axpy(E_[k]*Sk,S_[k],v);
      }
    }
  }

  void orthonormalize (MPS_type& c) const
  {
    for(size_t k = 0; k < N_-1; ++k) {
      double Sk = btas::Dotc(S_[k],c);
      btas::Axpy(-Sk,S_[k],c);
    }
    btas::Normalize(c);
  }

  const HOp_type* pLOp_;

  const HOp_type* pROp_;

  size_t N_;

  std::vector<double> E_;

  std::vector<MPS_type> H_;

  std::vector<MPS_type> S_;

};

//
// Davidson's precondition
//

template<size_t N, class Q>
void precondition (
  const double& eval,
  const btas::QSTArray<double,N,Q>& diag,
        btas::QSTArray<double,N,Q>& errv)
{
  for(auto ir = errv.begin(); ir != errv.end(); ++ir) {
    auto id = diag.find(ir->first);
    if(id != diag.end()) {
      auto irx = ir->second->begin();
      auto idx = id->second->begin();
      for(; irx != ir->second->end(); ++irx, ++idx) {
        double denm = eval - *idx;
        if(fabs(denm) < 1.0e-8) denm = 1.0e-8;
        *irx /= denm;
      }
    }
    else {
      btas::Dscal(1.0/eval,(*ir->second));
    }
  }
}

//
// Davidson eigen solver
//

template<size_t N, class Q>
double diagonalize (
  const Functor<N,Q>& multH,
  const btas::QSTArray<double,N,Q>& diag,
        btas::QSTArray<double,N,Q>& wfnc,
  const double& tolerance = 1.0e-8)
{
  Communicator world;

  int max_ritz = 10;
  double eval = 0.0;

  // reserve working space
  std::vector<btas::QSTArray<double,N,Q>> trial(max_ritz);
  std::vector<btas::QSTArray<double,N,Q>> sigma(max_ritz);

  btas::Copy(wfnc,trial[0]);

  int niter = 0;
  int iconv = 0;
  while(iconv < 1 && niter < 4) {
    pout << "\t\t\tmacro iteration [ " << std::setw(2) << niter << " ] " << std::endl;
    pout << "\t\t\t--------------------------------------------------" << std::endl;
    // to keep numerical stability
    btas::Normalize(trial[0]);
    sigma[0].clear();
    multH(trial[0],sigma[0]);
    for(int m = 1; m <= max_ritz; ++m) {
      // compute small Hamiltonian matrix
      btas::TArray<double,2> heff(m,m);
      btas::TArray<double,2> ovlp(m,m);
      for(int i = 0; i < m; ++i) {
        double hii;
#ifndef _SERIAL
        double hiitmp = btas::Dotc(trial[i],sigma[i]);
        boost::mpi::all_reduce(world,hiitmp,hii,std::plus<double>());
#else
        hii = btas::Dotc(trial[i],sigma[i]);
#endif
        heff(i,i) = hii;
        ovlp(i,i) = btas::Dotc(trial[i],trial[i]);
        for(int j = 0; j < i; ++j) {
          double hij;
#ifndef _SERIAL
          double hijtmp = btas::Dotc(trial[i],sigma[j]);
          boost::mpi::all_reduce(world,hijtmp,hij,std::plus<double>());
#else
          hij = btas::Dotc(trial[i],sigma[j]);
#endif
          heff(i,j) = hij;
          heff(j,i) = hij;
          double sij = btas::Dotc(trial[i],trial[j]);
          ovlp(i,j) = sij;
          ovlp(j,i) = sij;
        }
      }
      // solve eigenvalue problem to obtain Ritz value & vector
      btas::TArray<double,2> rvec;
      btas::TArray<double,1> rval;
      if(world.rank() == 0) {
        btas::Syev('V','U',heff,rval,rvec);
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,rval,0);
      boost::mpi::broadcast(world,rvec,0);
#endif
      eval = rval(0);
      pout << "\t\t\tmicro iteration [ " << std::setw(2) << m << " ] :: " << std::setprecision(16) << std::setw(24) << std::fixed << eval << std::endl;
      // rotate trial & sigma vectors by Ritz vector
      std::vector<btas::QSTArray<double,N,Q>> trial_save(m);
      std::vector<btas::QSTArray<double,N,Q>> sigma_save(m);
      for(int i = 0; i < m; ++i) {
        btas::Copy(trial[i],trial_save[i]);
        btas::Copy(sigma[i],sigma_save[i]);
        btas::Scal(rvec(i,i),trial[i]);
        btas::Scal(rvec(i,i),sigma[i]);
      }
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          if(i != j) {
            btas::Axpy(rvec(i,j),trial_save[i],trial[j]);
            btas::Axpy(rvec(i,j),sigma_save[i],sigma[j]);
          }
        }
      }
      // compute error vector
      btas::QSTArray<double,N,Q> errv;
      double rnorm;
#ifndef _SERIAL
      if(world.rank() != 0) {
        for(size_t i = 1; i < world.size(); ++i) {
          if(world.rank() == i) {
            world.send(0,i,sigma[0]);
          }
        }
      }
      else {
#endif
        btas::Copy(sigma[0],errv);
#ifndef _SERIAL
        for(size_t i = 1; i < world.size(); ++i) {
          btas::QSTArray<double,N,Q> temp;
          world.recv(i,i,temp);
          btas::Axpy(1.0,temp,errv);
        }
#endif
        btas::Axpy(-eval,trial[0],errv);
        rnorm = btas::Dotc(errv,errv);
#ifndef _SERIAL
      }
      boost::mpi::broadcast(world,rnorm,0);
#endif
      if(rnorm < tolerance) { ++iconv; break; }
      // solve correction equation
      if(m < max_ritz) {
        if(world.rank() == 0) {
          precondition(eval,diag,errv);
          for(int i = 0; i < m; ++i) {
            btas::Normalize(errv);
            btas::Orthogonalize(trial[i],errv);
          }
          btas::Normalize(errv);
          btas::Copy(errv,trial[m]);
        }
#ifndef _SERIAL
        boost::mpi::broadcast(world,trial[m],0);
#endif
        sigma[m].clear();
        multH(trial[m],sigma[m]);
      }
    }
    ++niter;
    pout << "\t\t\t--------------------------------------------------" << std::endl;
  }
  btas::Copy(trial[0],wfnc);

  multH.orthonormalize(wfnc);

  return eval;
}

} // namespace Davidson
} // namespace mpsxx

#endif // __MPSXX_DAVIDSON_DIAGONALIZE_HPP
