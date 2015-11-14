#ifndef __MPSXX_OPTIMIZE_HPP
#define __MPSXX_OPTIMIZE_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <btas/QSPARSE/QSTArray.h>
#include <time_stamp.h>

#include "mpidefs.h"

#include "renormalize.hpp"
#include "canonicalize.hpp"
#include "diagonalize.hpp"
#include "compute_diagonal_elements.hpp"

namespace mpsxx {

/// State-Specific DMRG optimization with 1-site algorithm
/// 1 forward or 1 backward sweep is performed.
/// \param forward perform forward sweep if true, and vice versa.
/// \param E energies up to the current root.
/// \param mpo a set of MPO tensors at this site, each of which is MPO in different group.
/// Normally, there is only 1 group per procs.
/// \param lHopr left-renormalized operators for each group.
/// \param rHopr right-renormalized operators for each group.
/// \param lSopr left-renormalized overlap.
/// \param rSopr right-renormalized overlap.
/// \param kMps left- or right-rotation matrix depending on \c forward .
/// \param lWfnc wavefunction at the left site (if fwd., this is system, else, environment).
/// \param rWfnc wavefunction at the right site (if fwd., this is environment, else, system).
/// \param tole tolerance in Davidson eigensolver.
/// \param M max. number of states to be kept in truncation (if 0, all (thin) singular vectors are kept).
/// \param noise random noise to be added in truncation.
template<class Q>
double optimize_onesite (
        bool forward,
  const std::vector<double>& E,
  const std::vector<btas::QSTArray<double,4,Q>>& mpo,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& lHopr,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& rHopr,
        std::vector<btas::QSTArray<double,2,Q>>& lSopr,
        std::vector<btas::QSTArray<double,2,Q>>& rSopr,
        std::vector<btas::QSTArray<double,3,Q>>& kMps,
        std::vector<btas::QSTArray<double,3,Q>>& lWfnc,
        std::vector<btas::QSTArray<double,3,Q>>& rWfnc,
  const double& tole, const int& M, const double& noise = 0.0)
{
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  Communicator world;

  size_t iroot = lWfnc.size()-1;

  time_stamp ts;

  double energy;

  if(forward) {
    pout << "\t\tcomputing diagonal elements..." << endl;
    btas::QSTArray<double,3,Q> diag(lWfnc[iroot].q(),lWfnc[iroot].qshape());
    compute_diagonal_elements(mpo,lHopr[iroot],rHopr[iroot],diag);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
    Davidson::Functor<3,Q> Fmult(E,mpo,lHopr,rHopr,lSopr,rSopr,lWfnc);
    energy = Davidson::diagonalize(Fmult,diag,lWfnc[iroot],tole);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    if(world.rank() == 0) {
      btas::QSTArray<double,3,Q> temp(lWfnc[iroot]);
      canonicalize(1,temp,lWfnc[iroot],rWfnc[iroot],M,noise);
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,lWfnc[iroot],0);
    boost::mpi::broadcast(world,rWfnc[iroot],0);
#endif
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = lWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      // renormalize for each operator group
      for(size_t g = 0; g < mpo.size(); ++g) {
        btas::QSTArray<double,3,Q> lHoprTmp;
        renormalize(1,mpo[g],lHopr[k][g],lWfnc[iroot],kMps[k],lHoprTmp);
        lHopr[k][g] = lHoprTmp;
      }

      btas::QSTArray<double,2,Q> lSoprTmp;
      renormalize(1,    lSopr[k],lWfnc[iroot],kMps[k],lSoprTmp);
      lSopr[k] = lSoprTmp;
    }
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    pout << "\t\tcomputing diagonal elements..." << endl;
    btas::QSTArray<double,3,Q> diag(rWfnc[iroot].q(),rWfnc[iroot].qshape());
    compute_diagonal_elements(mpo,lHopr[iroot],rHopr[iroot],diag);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
    Davidson::Functor<3,Q> Fmult(E,mpo,lHopr,rHopr,lSopr,rSopr,rWfnc);
    energy = Davidson::diagonalize(Fmult,diag,rWfnc[iroot],tole);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    if(world.rank() == 0) {
      btas::QSTArray<double,3,Q> temp(rWfnc[iroot]);
      canonicalize(0,temp,lWfnc[iroot],rWfnc[iroot],M,noise);
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,lWfnc[iroot],0);
    boost::mpi::broadcast(world,rWfnc[iroot],0);
#endif
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = rWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      // renormalize for each operator group
      for(size_t g = 0; g < mpo.size(); ++g) {
        btas::QSTArray<double,3,Q> rHoprTmp;
        renormalize(0,mpo[g],rHopr[k][g],rWfnc[iroot],kMps[k],rHoprTmp);
        rHopr[k][g] = rHoprTmp;
      }

      btas::QSTArray<double,2,Q> rSoprTmp;
      renormalize(0,    rSopr[k],rWfnc[iroot],kMps[k],rSoprTmp);
      rSopr[k] = rSoprTmp;
    }
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  pout << "\t\t--------------------------------------------------------------------------------" << endl;
  pout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

/// State-Specific DMRG optimization with 2-site algorithm
/// 1 forward or 1 backward sweep is performed.
/// \param forward perform forward sweep if true, and vice versa.
/// \param E energies up to the current root.
/// \param lmpo a set of MPO tensors at the left site, each of which is MPO in different group.
/// \param rmpo a set of MPO tensors at the right site, each of which is MPO in different group.
/// \param lHopr left-renormalized operators for each group.
/// \param rHopr right-renormalized operators for each group.
/// \param lSopr left-renormalized overlap.
/// \param rSopr right-renormalized overlap.
/// \param kMps left- or right-rotation matrix depending on \c forward .
/// \param lWfnc wavefunction at the left site (if fwd., this is system, else, environment).
/// \param rWfnc wavefunction at the right site (if fwd., this is environment, else, system).
/// \param tole tolerance in Davidson eigensolver.
/// \param M max. number of states to be kept in truncation (if 0, all (thin) singular vectors are kept).
/// \param noise random noise to be added in truncation.
template<class Q>
double optimize_twosite (
        bool forward,
  const std::vector<double>& E,
  const std::vector<btas::QSTArray<double,4,Q>>& lmpo,
  const std::vector<btas::QSTArray<double,4,Q>>& rmpo,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& lHopr,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& rHopr,
        std::vector<btas::QSTArray<double,2,Q>>& lSopr,
        std::vector<btas::QSTArray<double,2,Q>>& rSopr,
        std::vector<btas::QSTArray<double,3,Q>>& kMps,
        std::vector<btas::QSTArray<double,3,Q>>& lWfnc,
        std::vector<btas::QSTArray<double,3,Q>>& rWfnc,
  const double& tole, const int& M, const double& noise = 0.0)
{
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  Communicator world;

  size_t iroot = lWfnc.size()-1;

  time_stamp ts;

  pout << "\t\tcomputing 2-site wavefunction..." << endl;
  std::vector<btas::QSTArray<double,4,Q>> tWfnc(iroot+1);
  for(size_t k = 0; k <= iroot; ++k) {
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lWfnc[k],rWfnc[k],1.0,tWfnc[k]);
  }
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  pout << "\t\tcomputing diagonal elements..." << endl;
  btas::QSTArray<double,4,Q> diag(tWfnc[iroot].q(),tWfnc[iroot].qshape());
  compute_diagonal_elements(lmpo,rmpo,lHopr[iroot],rHopr[iroot],diag);
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  pout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
  Davidson::Functor<4,Q> Fmult(E,lmpo,rmpo,lHopr,rHopr,lSopr,rSopr,tWfnc);
  double energy = Davidson::diagonalize(Fmult,diag,tWfnc[iroot],tole);
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  if(forward) {
    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    if(world.rank() == 0) {
      canonicalize(1,tWfnc[iroot],lWfnc[iroot],rWfnc[iroot],M,noise);
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,lWfnc[iroot],0);
    boost::mpi::broadcast(world,rWfnc[iroot],0);
#endif
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = lWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      // renormalize for each operator group
      for(size_t g = 0; g < lmpo.size(); ++g) {
        btas::QSTArray<double,3,Q> lHoprTmp;
        renormalize(1,lmpo[g],lHopr[k][g],lWfnc[iroot],kMps[k],lHoprTmp);
        lHopr[k][g] = lHoprTmp;
      }

      btas::QSTArray<double,2,Q> lSoprTmp;
      renormalize(1,     lSopr[k],lWfnc[iroot],kMps[k],lSoprTmp);
      lSopr[k] = lSoprTmp;
    }
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    if(world.rank() == 0) {
      canonicalize(0,tWfnc[iroot],lWfnc[iroot],rWfnc[iroot],M,noise);
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,lWfnc[iroot],0);
    boost::mpi::broadcast(world,rWfnc[iroot],0);
#endif
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = rWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      // renormalize for each operator group
      for(size_t g = 0; g < rmpo.size(); ++g) {
        btas::QSTArray<double,3,Q> rHoprTmp;
        renormalize(0,rmpo[g],rHopr[k][g],rWfnc[iroot],kMps[k],rHoprTmp);
        rHopr[k][g] = rHoprTmp;
      }

      btas::QSTArray<double,2,Q> rSoprTmp;
      renormalize(0,     rSopr[k],rWfnc[iroot],kMps[k],rSoprTmp);
      rSopr[k] = rSoprTmp;
    }
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  pout << "\t\t--------------------------------------------------------------------------------" << endl;
  pout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

//
// With dot merge
//

/// Compute merged block operator with dot
template<class Q>
void compute_merged_operator (
        bool forward,
  const std::vector<btas::QSTArray<double,4,Q>>& xmpo,
  const std::vector<std::vector<btas::QSTArray<double,3,Q>>>& xHopr,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& mHopr,
  const std::vector<btas::QSTArray<double,2,Q>>& xSopr,
        std::vector<btas::QSTArray<double,2,Q>>& mSopr)
{
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  using btas::shape;

  size_t iroot = xHopr.size()-1;

  mHopr.resize(iroot+1);
  mSopr.resize(iroot+1);

  if(forward) {
    // Identity matrix acts on left site
    btas::QSTArray<double,2,Q>
    I(Q::zero(),btas::make_array(xmpo[0].qshape(1),xmpo[0].qshape(2)), btas::make_array(xmpo[0].dshape(1),xmpo[0].dshape(2)),1.0);

    for(size_t k = 0; k <= iroot; ++k) {
      const btas::Qshapes<Q>& q_l_bra = xHopr[k][0].qshape(0);
      const btas::Qshapes<Q>& q_n_bra = xmpo[0].qshape(1);
      const btas::Dshapes&    d_l_bra = xHopr[k][0].dshape(0);
      const btas::Dshapes&    d_n_bra = xmpo[0].dshape(1);

      const btas::Qshapes<Q>& q_l_ket = xHopr[k][0].qshape(2);
      const btas::Qshapes<Q>& q_n_ket = xmpo[0].qshape(2);
      const btas::Dshapes&    d_l_ket = xHopr[k][0].dshape(2);
      const btas::Dshapes&    d_n_ket = xmpo[0].dshape(2);

      btas::QSTmergeInfo<2,Q> M_q_bra(btas::make_array(q_l_bra,q_n_bra),btas::make_array(d_l_bra,d_n_bra));
      btas::QSTmergeInfo<2,Q> M_q_ket(btas::make_array(q_l_ket,q_n_ket),btas::make_array(d_l_ket,d_n_ket));

      mHopr[k].resize(xmpo.size());

      // merge operators for each group
      for(size_t g = 0; g < xmpo.size(); ++g) {
        btas::QSTArray<double,5,Q> lHopr;
        btas::Contract(1.0,xHopr[k][g],shape(0,1,2),xmpo[g],shape(1,3,4,5),1.0,lHopr,shape(0,3,5,2,4));

        btas::QSTArray<double,4,Q> lHscr;
        btas::QSTmerge(M_q_bra,lHopr,lHscr);
        btas::QSTmerge(lHscr,M_q_ket,mHopr[k][g]);
      }

      // merge overlap matrix
      btas::QSTArray<double,4,Q> lSscr;
      btas::Contract(1.0,xSopr[k],shape(0,1),I,shape(2,3),1.0,lSscr,shape(0,2,1,3));
      btas::QSTmerge(M_q_bra,lSscr,M_q_ket,mSopr[k]);
    }
  }
  else {
    // Identity matrix acts on right site
    btas::QSTArray<double,2,Q>
    I(Q::zero(),btas::make_array(xmpo[0].qshape(1),xmpo[0].qshape(2)), btas::make_array(xmpo[0].dshape(1),xmpo[0].dshape(2)),1.0);

    for(size_t k = 0; k <= iroot; ++k) {
      const btas::Qshapes<Q>& q_n_bra = xmpo[0].qshape(1);
      const btas::Qshapes<Q>& q_r_bra = xHopr[k][0].qshape(0);
      const btas::Dshapes&    d_n_bra = xmpo[0].dshape(1);
      const btas::Dshapes&    d_r_bra = xHopr[k][0].dshape(0);

      const btas::Qshapes<Q>& q_n_ket = xmpo[0].qshape(2);
      const btas::Qshapes<Q>& q_r_ket = xHopr[k][0].qshape(2);
      const btas::Dshapes&    d_n_ket = xmpo[0].dshape(1);
      const btas::Dshapes&    d_r_ket = xHopr[k][0].dshape(0);

      btas::QSTmergeInfo<2,Q> M_q_bra(btas::make_array(q_n_bra,q_r_bra),btas::make_array(d_n_bra,d_r_bra));
      btas::QSTmergeInfo<2,Q> M_q_ket(btas::make_array(q_n_ket,q_r_ket),btas::make_array(d_n_ket,d_r_ket));

      mHopr[k].resize(xmpo.size());

      // merge operators for each group
      for(size_t g = 0; g < xmpo.size(); ++g) {
        btas::QSTArray<double,5,Q> rHopr;
        btas::Contract(1.0,xmpo[g],shape(0,1,2,3),xHopr[k][g],shape(4,3,5),1.0,rHopr,shape(1,4,0,2,5));

        btas::QSTArray<double,4,Q> rHscr;
        btas::QSTmerge(M_q_bra,rHopr,rHscr);
        btas::QSTmerge(rHscr,M_q_ket,mHopr[k][g]);
      }

      // merge overlap matrix
      btas::QSTArray<double,4,Q> rSscr;
      btas::Contract(1.0,I,shape(0,1),xSopr[k],shape(2,3),1.0,rSscr,shape(0,2,1,3));
      btas::QSTmerge(M_q_bra,rSscr,M_q_ket,mSopr[k]);
    }
  }
}

/***
template<class Q>
double M_optimize_onesite (
        bool forward,
  const std::vector<double>& E,
  const btas::QSTArray<double,4,Q>& mpo,
        std::vector<btas::QSTArray<double,3,Q>>& lHopr,
        std::vector<btas::QSTArray<double,3,Q>>& rHopr,
        std::vector<btas::QSTArray<double,2,Q>>& lSopr,
        std::vector<btas::QSTArray<double,2,Q>>& rSopr,
        std::vector<btas::QSTArray<double,3,Q>>& kMps,
        std::vector<btas::QSTArray<double,3,Q>>& lWfnc,
        std::vector<btas::QSTArray<double,3,Q>>& rWfnc,
               const double& tole, const int& M = 0)
{
  using std::endl;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  size_t iroot = lWfnc.size()-1;

  std::vector<btas::QSTArray<double,3,Q>> M_lHopr(iroot+1);
  std::vector<btas::QSTArray<double,3,Q>> M_rHopr(iroot+1);
  std::vector<btas::QSTArray<double,2,Q>> M_lSopr(iroot+1);
  std::vector<btas::QSTArray<double,2,Q>> M_rSopr(iroot+1);
  std::vector<btas::QSTArray<double,2,Q>> M_tWfnc(iroot+1);

  std::vector<btas::QSTmergeInfo<2,Q>> M_q_ket(iroot+1);

  time_stamp ts;

  pout << "\t\tconstructing merged super blocks..." << endl;
  if(forward) {
    for(size_t k = 0; k <= iroot; ++k) {
      const btas::Qshapes<Q>& q_l_ket =-lHopr[k].qshape(2);
      const btas::Qshapes<Q>& q_n_ket =-mpo.qshape(2);
      const btas::Dshapes&    d_l_ket = lHopr[k].dshape(0);
      const btas::Dshapes&    d_n_ket = mpo.dshape(1);
      M_q_ket.reset(btas::make_array(q_l_ket,q_n_ket),btas::make_array(d_l_ket,d_n_ket));

      btas::QSTmerge(M_q_ket,lWfnc,M_tWfnc);

      make_merge(1,mpo,lHopr,M_lHopr);
      M_rHopr.reference(rHopr);
    }
  }
  else {
    const btas::Qshapes<Q>& q_n_ket =-mpo.qshape(2);
    const btas::Qshapes<Q>& q_r_ket =-rHopr.qshape(2);
    const btas::Dshapes&    d_n_ket = mpo.dshape(1);
    const btas::Dshapes&    d_r_ket = rHopr.dshape(0);
    M_q_ket.reset(btas::make_array(q_n_ket,q_r_ket),btas::make_array(d_n_ket,d_r_ket));

    btas::QSTmerge(lWfnc,M_q_ket,M_tWfnc);

    M_lHopr.reference(lHopr);
    make_merge(0,mpo,rHopr,M_rHopr);
  }
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  pout << "\t\tcomputing diagonal elements..." << endl;
  btas::QSTArray<double,2,Q> diag(M_tWfnc.q(),M_tWfnc.qshape());
  compute_diagonal_elements(M_lHopr,M_rHopr,diag);
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  pout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
  Davidson::Functor<2,Q> Fmult = boost::bind(compute_sigmavector<Q>,M_lHopr,M_rHopr,_1,_2);
  double energy = Davidson::diagonalize(Fmult,diag,M_tWfnc,tole);
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  if(forward) {
    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    btas::QSTArray<double,2,Q> lmps_mg;
    btas::QSTArray<double,2,Q> gmat;
    canonicalize(1,M_tWfnc,lmps_mg,gmat,M,noise);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(M_q_ket,lmps_mg,lWfnc);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\tcomputing guess wavefunction to the next..." << endl;
    btas::QSTArray<double,3,Q> rWfncTmp;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,gmat,rWfnc,1.0,rWfncTmp);
    rWfnc = rWfncTmp;
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    lHopr.clear();
    renormalize(1,M_lHopr,lmps_mg,lmps_mg,lHopr);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    btas::QSTArray<double,2,Q> rmps_mg;
    btas::QSTArray<double,2,Q> gmat;
    canonicalize(0,M_tWfnc,rmps_mg,gmat,M,noise);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(rmps_mg,M_q_ket,lWfnc);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\tcomputing guess wavefunction to the next..." << endl;
    btas::QSTArray<double,3,Q> rWfncTmp;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,rWfnc,gmat,1.0,rWfncTmp);
    rWfnc = rWfncTmp;
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    rHopr.clear();
    renormalize(0,M_rHopr,rmps_mg,rmps_mg,rHopr);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  pout << "\t\t--------------------------------------------------------------------------------" << endl;
  pout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}
***/

template<class Q>
double optimize_twosite_merged (
        bool forward,
  const std::vector<double>& E,
  const std::vector<btas::QSTArray<double,4,Q>>& lmpo,
  const std::vector<btas::QSTArray<double,4,Q>>& rmpo,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& lHopr,
        std::vector<std::vector<btas::QSTArray<double,3,Q>>>& rHopr,
        std::vector<btas::QSTArray<double,2,Q>>& lSopr,
        std::vector<btas::QSTArray<double,2,Q>>& rSopr,
        std::vector<btas::QSTArray<double,3,Q>>& kMps,
        std::vector<btas::QSTArray<double,3,Q>>& lWfnc,
        std::vector<btas::QSTArray<double,3,Q>>& rWfnc,
  const double& tole, const int& M, const double& noise = 0.0)
{
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  Communicator world;

  size_t iroot = lWfnc.size()-1;

  time_stamp ts;

  pout << "\t\tconstructing merged super blocks..." << endl;

  // merge left block operators
  std::vector<std::vector<btas::QSTArray<double,3,Q>>> M_lHopr;
  std::vector<btas::QSTArray<double,2,Q>> M_lSopr;

  compute_merged_operator(1,lmpo,lHopr,M_lHopr,lSopr,M_lSopr);

  // merge right block operators
  std::vector<std::vector<btas::QSTArray<double,3,Q>>> M_rHopr;
  std::vector<btas::QSTArray<double,2,Q>> M_rSopr;

  compute_merged_operator(0,rmpo,rHopr,M_rHopr,rSopr,M_rSopr);

  // merged wavefunction
  std::vector<btas::QSTArray<double,2,Q>> M_tWfnc(iroot+1);
  std::vector<btas::QSTArray<double,2,Q>> M_kMps(iroot+1);

  btas::QSTmergeInfo<2,Q> q_lM_ket;
  btas::QSTmergeInfo<2,Q> q_rM_ket;

  for(size_t k = 0; k <= iroot; ++k) {
    const btas::Qshapes<Q>& q_l_ket =-lHopr[k][0].qshape(2);
    const btas::Qshapes<Q>& q_m_ket =-lmpo[0].qshape(2);
    const btas::Dshapes&    d_l_ket = lHopr[k][0].dshape(0);
    const btas::Dshapes&    d_m_ket = lmpo[0].dshape(1);
    q_lM_ket.reset(btas::make_array(q_l_ket,q_m_ket),btas::make_array(d_l_ket,d_m_ket));

    const btas::Qshapes<Q>& q_n_ket =-rmpo[0].qshape(2);
    const btas::Qshapes<Q>& q_r_ket =-rHopr[k][0].qshape(2);
    const btas::Dshapes&    d_n_ket = rmpo[0].dshape(1);
    const btas::Dshapes&    d_r_ket = rHopr[k][0].dshape(0);
    q_rM_ket.reset(btas::make_array(q_n_ket,q_r_ket),btas::make_array(d_n_ket,d_r_ket));

    btas::QSTArray<double,4,Q> tWfnc;
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lWfnc[k],rWfnc[k],1.0,tWfnc);
    btas::QSTmerge(q_lM_ket,tWfnc,q_rM_ket,M_tWfnc[k]);

    if(k < iroot) {
      if(forward)
        btas::QSTmerge(q_lM_ket,kMps[k],M_kMps[k]);
      else
        btas::QSTmerge(kMps[k],q_rM_ket,M_kMps[k]);
    }
  }
  // NOTE: on exit the loop, q_lM_ket and q_rM_ket are stored for iroot to expand the wavefunction.

  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  pout << "\t\tcomputing diagonal elements..." << endl;
  btas::QSTArray<double,2,Q> diag(M_tWfnc[iroot].q(),M_tWfnc[iroot].qshape());
  compute_diagonal_elements(M_lHopr[iroot],M_rHopr[iroot],diag);
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  pout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
  Davidson::Functor<2,Q> Fmult(E,M_lHopr,M_rHopr,M_lSopr,M_rSopr,M_tWfnc);
  double energy = Davidson::diagonalize(Fmult,diag,M_tWfnc[iroot],tole);
  pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  btas::QSTArray<double,2,Q> M_lWfnc;
  btas::QSTArray<double,2,Q> M_rWfnc;

  if(forward) {
    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    if(world.rank() == 0) {
      canonicalize(1,M_tWfnc[iroot],M_lWfnc,M_rWfnc,M,noise);
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,M_lWfnc,0);
    boost::mpi::broadcast(world,M_rWfnc,0);
#endif
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    M_kMps[iroot] = M_lWfnc;
    for(size_t k = 0; k <= iroot; ++k) {
      // renormalize for each operator group
      for(size_t g = 0; g < lmpo.size(); ++g) {
        btas::QSTArray<double,3,Q> lHoprTmp;
        renormalize(1,M_lHopr[k][g],M_lWfnc,M_kMps[k],lHoprTmp);
        lHopr[k][g] = lHoprTmp;
      }

      btas::QSTArray<double,2,Q> lSoprTmp;
      renormalize(1,M_lSopr[k],M_lWfnc,M_kMps[k],lSoprTmp);
      lSopr[k] = lSoprTmp;
    }
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(q_lM_ket,M_lWfnc,lWfnc[iroot]);
    btas::QSTexpand(M_rWfnc,q_rM_ket,rWfnc[iroot]);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    pout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    if(world.rank() == 0) {
      canonicalize(0,M_tWfnc[iroot],M_lWfnc,M_rWfnc,M,noise);
    }
#ifndef _SERIAL
    boost::mpi::broadcast(world,M_lWfnc,0);
    boost::mpi::broadcast(world,M_rWfnc,0);
#endif
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\trenormalizing operators to the next..." << endl;
    M_kMps[iroot] = M_rWfnc;
    for(size_t k = 0; k <= iroot; ++k) {
      // renormalize for each operator group
      for(size_t g = 0; g < rmpo.size(); ++g) {
        btas::QSTArray<double,3,Q> rHoprTmp;
        renormalize(0,M_rHopr[k][g],M_rWfnc,M_kMps[k],rHoprTmp);
        rHopr[k][g] = rHoprTmp;
      }

      btas::QSTArray<double,2,Q> rSoprTmp;
      renormalize(0,M_rSopr[k],M_rWfnc,M_kMps[k],rSoprTmp);
      rSopr[k] = rSoprTmp;
    }
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    pout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(q_lM_ket,M_lWfnc,lWfnc[iroot]);
    btas::QSTexpand(M_rWfnc,q_rM_ket,rWfnc[iroot]);
    pout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  pout << "\t\t--------------------------------------------------------------------------------" << endl;
  pout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

} // namespace mpsxx

#endif // __MPSXX_OPTIMIZE_HPP
