#ifndef __MPSXX_OPTIMIZE_HPP
#define __MPSXX_OPTIMIZE_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <legacy/QSPARSE/QSTArray.h>
#include <time_stamp.h>

#include "renormalize.hpp"
#include "canonicalize.hpp"
#include "diagonalize.hpp"
#include "compute_diagonal_elements.hpp"

namespace mpsxx {

template<class Q>
double optimize_onesite (
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
  const double& T, const int& M)
{
  using std::cout;
  using std::endl;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  size_t iroot = lWfnc.size()-1;

  time_stamp ts;

  double energy;

  if(forward) {
    cout << "\t\tcomputing diagonal elements..." << endl;
    btas::QSTArray<double,3,Q> diag(lWfnc[iroot].q(),lWfnc[iroot].qshape());
    compute_diagonal_elements(mpo,lHopr[iroot],rHopr[iroot],diag);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
    Davidson::Functor<3,Q> Fmult(E,mpo,lHopr,rHopr,lSopr,rSopr,lWfnc);
    energy = Davidson::diagonalize(Fmult,diag,lWfnc[iroot]);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    btas::QSTArray<double,3,Q> temp(lWfnc[iroot]);
    canonicalize(1,temp,lWfnc[iroot],rWfnc[iroot],M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = lWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      btas::QSTArray<double,3,Q> lHoprTmp;
      renormalize(1,mpo,lHopr[k],lWfnc[iroot],kMps[k],lHoprTmp);
      lHopr[k] = lHoprTmp;

      btas::QSTArray<double,2,Q> lSoprTmp;
      renormalize(1,    lSopr[k],lWfnc[iroot],kMps[k],lSoprTmp);
      lSopr[k] = lSoprTmp;
    }
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tcomputing diagonal elements..." << endl;
    btas::QSTArray<double,3,Q> diag(rWfnc[iroot].q(),rWfnc[iroot].qshape());
    compute_diagonal_elements(mpo,lHopr[iroot],rHopr[iroot],diag);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
    Davidson::Functor<3,Q> Fmult(E,mpo,lHopr,rHopr,lSopr,rSopr,rWfnc);
    energy = Davidson::diagonalize(Fmult,diag,rWfnc[iroot]);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    btas::QSTArray<double,3,Q> temp(rWfnc[iroot]);
    canonicalize(0,temp,lWfnc[iroot],rWfnc[iroot],M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = rWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      btas::QSTArray<double,3,Q> rHoprTmp;
      renormalize(0,mpo,rHopr[k],rWfnc[iroot],kMps[k],rHoprTmp);
      rHopr[k] = rHoprTmp;

      btas::QSTArray<double,2,Q> rSoprTmp;
      renormalize(0,    rSopr[k],rWfnc[iroot],kMps[k],rSoprTmp);
      rSopr[k] = rSoprTmp;
    }
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

template<class Q>
double optimize_twosite (
        bool forward,
  const std::vector<double>& E,
  const btas::QSTArray<double,4,Q>& lmpo,
  const btas::QSTArray<double,4,Q>& rmpo,
        std::vector<btas::QSTArray<double,3,Q>>& lHopr,
        std::vector<btas::QSTArray<double,3,Q>>& rHopr,
        std::vector<btas::QSTArray<double,2,Q>>& lSopr,
        std::vector<btas::QSTArray<double,2,Q>>& rSopr,
        std::vector<btas::QSTArray<double,3,Q>>& kMps,
        std::vector<btas::QSTArray<double,3,Q>>& lWfnc,
        std::vector<btas::QSTArray<double,3,Q>>& rWfnc,
  const double& T, const int& M)
{
  using std::cout;
  using std::endl;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  size_t iroot = lWfnc.size()-1;

  time_stamp ts;

  cout << "\t\tcomputing 2-site wavefunction..." << endl;
  std::vector<btas::QSTArray<double,4,Q>> tWfnc(iroot+1);
  for(size_t k = 0; k <= iroot; ++k) {
    btas::Gemm(btas::NoTrans,btas::NoTrans,1.0,lWfnc[k],rWfnc[k],1.0,tWfnc[k]);
  }
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\tcomputing diagonal elements..." << endl;
  btas::QSTArray<double,4,Q> diag(tWfnc[iroot].q(),tWfnc[iroot].qshape());
  compute_diagonal_elements(lmpo,rmpo,lHopr[iroot],rHopr[iroot],diag);
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
  Davidson::Functor<4,Q> Fmult(E,lmpo,rmpo,lHopr,rHopr,lSopr,rSopr,tWfnc);
  double energy = Davidson::diagonalize(Fmult,diag,tWfnc[iroot]);
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    canonicalize(1,tWfnc[iroot],lWfnc[iroot],rWfnc[iroot],M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = lWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      btas::QSTArray<double,3,Q> lHoprTmp;
      renormalize(1,lmpo,lHopr[k],lWfnc[iroot],kMps[k],lHoprTmp);
      lHopr[k] = lHoprTmp;

      btas::QSTArray<double,2,Q> lSoprTmp;
      renormalize(1,     lSopr[k],lWfnc[iroot],kMps[k],lSoprTmp);
      lSopr[k] = lSoprTmp;
    }
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    canonicalize(0,tWfnc[iroot],lWfnc[iroot],rWfnc[iroot],M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    kMps[iroot] = rWfnc[iroot]; // FIXME: this is extra operation.
    for(size_t k = 0; k <= iroot; ++k) {
      btas::QSTArray<double,3,Q> rHoprTmp;
      renormalize(0,rmpo,rHopr[k],rWfnc[iroot],kMps[k],rHoprTmp);
      rHopr[k] = rHoprTmp;

      btas::QSTArray<double,2,Q> rSoprTmp;
      renormalize(0,     rSopr[k],rWfnc[iroot],kMps[k],rSoprTmp);
      rSopr[k] = rSoprTmp;
    }
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

//
// With dot merge
//

/// Compute merged block with dot
template<class Q>
void compute_merged_block (
        bool forward,
  const btas::QSTArray<double,4,Q>& mpo0,
        btas::QSTArray<double,3,Q>& opr0,
        btas::QSTArray<double,3,Q>& oprm)
{
  using std::cout;
  using std::endl;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  using btas::shape;

  if(forward) {
    const btas::Qshapes<Q>& q_l_bra = opr0.qshape(0);
    const btas::Qshapes<Q>& q_n_bra = mpo0.qshape(1);
    const btas::Dshapes&    d_l_bra = opr0.dshape(0);
    const btas::Dshapes&    d_n_bra = mpo0.dshape(1);

    const btas::Qshapes<Q>& q_l_ket = opr0.qshape(2);
    const btas::Qshapes<Q>& q_n_ket = mpo0.qshape(2);
    const btas::Dshapes&    d_l_ket = opr0.dshape(0);
    const btas::Dshapes&    d_n_ket = mpo0.dshape(1);

    btas::QSTArray<double,5,Q> lopr;
    btas::Contract(1.0,opr0,shape(0,1,2),mpo0,shape(1,3,4,5),1.0,lopr,shape(0,3,5,2,4));

    btas::QSTmergeInfo<2,Q> q_mg_bra(btas::make_array(q_l_bra,q_n_bra),btas::make_array(d_l_bra,d_n_bra));
    btas::QSTmergeInfo<2,Q> q_mg_ket(btas::make_array(q_l_ket,q_n_ket),btas::make_array(d_l_ket,d_n_ket));

    btas::QSTArray<double,4,Q> oprx;
    btas::QSTmerge(q_mg_bra,lopr,oprx);
    btas::QSTmerge(oprx,q_mg_ket,oprm);
  }
  else {
    const btas::Qshapes<Q>& q_n_bra = mpo0.qshape(1);
    const btas::Qshapes<Q>& q_r_bra = opr0.qshape(0);
    const btas::Dshapes&    d_n_bra = mpo0.dshape(1);
    const btas::Dshapes&    d_r_bra = opr0.dshape(0);

    const btas::Qshapes<Q>& q_n_ket = mpo0.qshape(2);
    const btas::Qshapes<Q>& q_r_ket = opr0.qshape(2);
    const btas::Dshapes&    d_n_ket = mpo0.dshape(1);
    const btas::Dshapes&    d_r_ket = opr0.dshape(0);

    btas::QSTArray<double,5,Q> ropr;
    btas::Contract(1.0,mpo0,shape(0,1,2,3),opr0,shape(4,3,5),1.0,ropr,shape(1,4,0,2,5));

    btas::QSTmergeInfo<2,Q> q_mg_bra(btas::make_array(q_n_bra,q_r_bra),btas::make_array(d_n_bra,d_r_bra));
    btas::QSTmergeInfo<2,Q> q_mg_ket(btas::make_array(q_n_ket,q_r_ket),btas::make_array(d_n_ket,d_r_ket));

    btas::QSTArray<double,4,Q> oprx;
    btas::QSTmerge(q_mg_bra,ropr,oprx);
    btas::QSTmerge(oprx,q_mg_ket,oprm);
  }
}

/***
template<class Q>
double optimize_onesite_merged
(bool forward, const btas::QSTArray<double,4,Q>& mpo0,
                     btas::QSTArray<double,3,Q>& lopr,
                     btas::QSTArray<double,3,Q>& ropr,
                     btas::QSTArray<double,3,Q>& lWfnc,
                     btas::QSTArray<double,3,Q>& rWfnc,
               const double& T, const int& M = 0)
{
  using std::cout;
  using std::endl;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  btas::QSTArray<double,3,Q> lopr_mg;
  btas::QSTArray<double,3,Q> ropr_mg;
  btas::QSTArray<double,2,Q> wfnc_mg;
  btas::QSTmergeInfo<2,Q> q_mg_ket;

  time_stamp ts;

  cout << "\t\tconstructing merged super blocks..." << endl;
  if(forward) {
    const btas::Qshapes<Q>& q_l_ket =-lopr.qshape(2);
    const btas::Qshapes<Q>& q_n_ket =-mpo0.qshape(2);
    const btas::Dshapes&    d_l_ket = lopr.dshape(0);
    const btas::Dshapes&    d_n_ket = mpo0.dshape(1);
    q_mg_ket.reset(btas::make_array(q_l_ket, q_n_ket), btas::make_array(d_l_ket, d_n_ket));

    btas::QSTmerge(q_mg_ket, lWfnc, wfnc_mg);

    compute_merged_block(1, mpo0, lopr, lopr_mg);
    ropr_mg.reference(ropr);
  }
  else {
    const btas::Qshapes<Q>& q_n_ket =-mpo0.qshape(2);
    const btas::Qshapes<Q>& q_r_ket =-ropr.qshape(2);
    const btas::Dshapes&    d_n_ket = mpo0.dshape(1);
    const btas::Dshapes&    d_r_ket = ropr.dshape(0);
    q_mg_ket.reset(btas::make_array(q_n_ket, q_r_ket), btas::make_array(d_n_ket, d_r_ket));

    btas::QSTmerge(lWfnc, q_mg_ket, wfnc_mg);

    lopr_mg.reference(lopr);
    compute_merged_block(0, mpo0, ropr, ropr_mg);
  }
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\tcomputing diagonal elements..." << endl;
  btas::QSTArray<double,2,Q> diag(wfnc_mg.q(), wfnc_mg.qshape());
  compute_diagonal_elements(lopr_mg, ropr_mg, diag);
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
  Davidson::Functor<2,Q> Fmult = boost::bind(compute_sigmavector<Q>, lopr_mg, ropr_mg, _1, _2);
  double energy = Davidson::diagonalize(Fmult, diag, wfnc_mg);
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    btas::QSTArray<double,2,Q> lmps_mg;
    btas::QSTArray<double,2,Q> gmat;
    canonicalize(1, wfnc_mg, lmps_mg, gmat, M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(q_mg_ket, lmps_mg, lWfnc);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tcomputing guess wavefunction to the next..." << endl;
    btas::QSTArray<double,3,Q> rWfncTmp;
    btas::Gemm(btas::NoTrans, btas::NoTrans, 1.0, gmat, rWfnc, 1.0, rWfncTmp);
    rWfnc = rWfncTmp;
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    lopr.clear();
    renormalize(1, lopr_mg, lmps_mg, lmps_mg, lopr);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    btas::QSTArray<double,2,Q> rmps_mg;
    btas::QSTArray<double,2,Q> gmat;
    canonicalize(0, wfnc_mg, rmps_mg, gmat, M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(rmps_mg, q_mg_ket, lWfnc);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tcomputing guess wavefunction to the next..." << endl;
    btas::QSTArray<double,3,Q> rWfncTmp;
    btas::Gemm(btas::NoTrans, btas::NoTrans, 1.0, rWfnc, gmat, 1.0, rWfncTmp);
    rWfnc = rWfncTmp;
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    ropr.clear();
    renormalize(0, ropr_mg, rmps_mg, rmps_mg, ropr);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

template<class Q>
double optimize_twosite_merged (
        bool forward,
  const btas::QSTArray<double,4,Q>& lmpo,
               const btas::QSTArray<double,4,Q>& rmpo,
                     btas::QSTArray<double,3,Q>& lopr,
                     btas::QSTArray<double,3,Q>& ropr,
                     btas::QSTArray<double,3,Q>& lwfn,
                     btas::QSTArray<double,3,Q>& rwfn,
               const double& T, const int& M = 0)
{
  using std::cout;
  using std::endl;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  time_stamp ts;

  cout << "\t\tconstructing merged super blocks..." << endl;

  const btas::Qshapes<Q>& q_l_ket =-lopr.qshape(2);
  const btas::Qshapes<Q>& q_m_ket =-lmpo.qshape(2);
  const btas::Dshapes&    d_l_ket = lopr.dshape(0);
  const btas::Dshapes&    d_m_ket = lmpo.dshape(1);
  btas::QSTmergeInfo<2,Q> q_lmg_ket(btas::make_array(q_l_ket, q_m_ket), btas::make_array(d_l_ket, d_m_ket));

  btas::QSTArray<double,3,Q> lopr_mg;
  compute_merged_block(1, lmpo, lopr, lopr_mg);

  const btas::Qshapes<Q>& q_n_ket =-rmpo.qshape(2);
  const btas::Qshapes<Q>& q_r_ket =-ropr.qshape(2);
  const btas::Dshapes&    d_n_ket = rmpo.dshape(1);
  const btas::Dshapes&    d_r_ket = ropr.dshape(0);
  btas::QSTmergeInfo<2,Q> q_rmg_ket(btas::make_array(q_n_ket, q_r_ket), btas::make_array(d_n_ket, d_r_ket));

  btas::QSTArray<double,3,Q> ropr_mg;
  compute_merged_block(0, rmpo, ropr, ropr_mg);

  btas::QSTArray<double,2,Q> wfnc_mg;
  {
    btas::QSTArray<double,4,Q> tWfnc;
    btas::Gemm(btas::NoTrans, btas::NoTrans, 1.0, lwfn, rwfn, 1.0, tWfnc);
    btas::QSTmerge(q_lmg_ket, tWfnc, q_rmg_ket, wfnc_mg);
  }

  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\tcomputing diagonal elements..." << endl;
  btas::QSTArray<double,2,Q> diag(wfnc_mg.q(), wfnc_mg.qshape());
  compute_diagonal_elements(lopr_mg, ropr_mg, diag);
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << endl;
  Davidson::Functor<2,Q> Fmult = boost::bind(compute_sigmavector<Q>, lopr_mg, ropr_mg, _1, _2);
  double energy = Davidson::diagonalize(Fmult, diag, wfnc_mg);
  cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  btas::QSTArray<double,2,Q> lwfn_mg;
  btas::QSTArray<double,2,Q> rwfn_mg;

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    canonicalize(1, wfnc_mg, lwfn_mg, rwfn_mg, M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(q_lmg_ket, lwfn_mg, lwfn);
    btas::QSTexpand(rwfn_mg, q_rmg_ket, rwfn);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    lopr.clear();
    renormalize(1, lopr_mg, lwfn_mg, lwfn_mg, lopr);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << endl;
    canonicalize(0, wfnc_mg, rwfn_mg, lwfn_mg, M);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << endl;
    btas::QSTexpand(q_lmg_ket, lwfn_mg, lwfn);
    btas::QSTexpand(rwfn_mg, q_rmg_ket, rwfn);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << endl;
    ropr.clear();
    renormalize(0, ropr_mg, rwfn_mg, rwfn_mg, ropr);
    cout << "\t\tdone ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

***/

} // namespace mpsxx

#endif // __MPSXX_OPTIMIZE_HPP
