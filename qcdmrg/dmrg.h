#ifndef _MPSXX_CXX11_DMRG_H
#define _MPSXX_CXX11_DMRG_H 1

#include <iostream>
#include <iomanip>

#include <btas/QSPARSE/QSDArray.h>
#include <time_stamp.h>

#include <mpsxx.h>

#include <MpOperators.h>
#include <MpStates.h>

#include <driver/canonicalize.h>
#include <driver/renormalize.h>
#include <driver/guesswave.h>
#include <driver/diagelements.h>
#include <driver/sigmavector.h>
#include <driver/davidson.h>
#include <driver/fileio.h>

#include "input.h"

namespace mpsxx {

enum DMRGALGORITHM { ONESITE, TWOSITE };

template<class Q>
double optimize_onesite
(bool forward, const btas::QSDArray<4, Q>& mpo0,
                     btas::QSDArray<3, Q>& lopr,
                     btas::QSDArray<3, Q>& ropr,
                     btas::QSDArray<3, Q>& wfn0,
                     btas::QSDArray<3, Q>& wfn1,
               const double& T, const int& M = 0)
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  time_stamp ts;

  cout << "\t\tcomputing diagonal elements..." << flush;
  btas::QSDArray<3, Q> diag(wfn0.q(), wfn0.qshape());
  compute_diagonal_elements(mpo0, lopr, ropr, diag);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << flush;
  davidson::Functor<3, Q> f_sigmavector = boost::bind(compute_sigmavector<Q>, mpo0, lopr, ropr, _1, _2);
  double energy = davidson::diagonalize(f_sigmavector, diag, wfn0);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    btas::QSDArray<3, Q> lmps;
    btas::QSDArray<2, Q> gaug;
    canonicalize(1, wfn0, lmps, gaug, M);
    wfn0 = lmps;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tcompute guess wavefunction to the next..." << flush;
    btas::QSDArray<3, Q> wfn1_tmp;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, gaug, wfn1, 1.0, wfn1_tmp);
    wfn1 = wfn1_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    btas::QSDArray<3, Q> lopr_tmp;
    renormalize(1, mpo0, lopr, lmps, lmps, lopr_tmp);
    lopr = lopr_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    btas::QSDArray<3, Q> rmps;
    btas::QSDArray<2, Q> gaug;
    canonicalize(0, wfn0, rmps, gaug, M);
    wfn0 = rmps;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tcompute guess wavefunction to the next..." << flush;
    btas::QSDArray<3, Q> wfn1_tmp;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, wfn1, gaug, 1.0, wfn1_tmp);
    wfn1 = wfn1_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    btas::QSDArray<3, Q> ropr_tmp;
    renormalize(0, mpo0, ropr, rmps, rmps, ropr_tmp);
    ropr = ropr_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

template<class Q>
double optimize_twosite
(bool forward, const btas::QSDArray<4, Q>& lmpo,
               const btas::QSDArray<4, Q>& rmpo,
                     btas::QSDArray<3, Q>& lopr,
                     btas::QSDArray<3, Q>& ropr,
                     btas::QSDArray<3, Q>& lwfn,
                     btas::QSDArray<3, Q>& rwfn,
               const double& T, const int& M = 0)
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  time_stamp ts;

  cout << "\t\tcomputing 2-site wavefunction..." << flush;
  btas::QSDArray<4, Q> wfn2;
  btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, lwfn, rwfn, 1.0, wfn2);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\tcomputing diagonal elements..." << flush;
  btas::QSDArray<4, Q> diag(wfn2.q(), wfn2.qshape());
  compute_diagonal_elements(lmpo, rmpo, lopr, ropr, diag);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << flush;
  davidson::Functor<4, Q> f_sigmavector;
  f_sigmavector = boost::bind(compute_sigmavector<Q>, lmpo, rmpo, lopr, ropr, _1, _2);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  double energy = davidson::diagonalize(f_sigmavector, diag, wfn2);

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    canonicalize(1, wfn2, lwfn, rwfn, M);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    btas::QSDArray<3, Q> lopr_tmp;
    renormalize(1, lmpo, lopr, lwfn, lwfn, lopr_tmp);
    lopr = lopr_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    canonicalize(0, wfn2, rwfn, lwfn, M);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    btas::QSDArray<3, Q> ropr_tmp;
    renormalize(0, rmpo, ropr, rwfn, rwfn, ropr_tmp);
    ropr = ropr_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

//
// Merged block version
//

template<class Q>
void compute_merged_block
(bool forward, const btas::QSDArray<4, Q>& mpo0,
                     btas::QSDArray<3, Q>& opr0,
                     btas::QSDArray<3, Q>& oprm)
{
  using std::cout;
  using std::endl;
  using std::flush;
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

    btas::QSDArray<5, Q> lopr;
    btas::QSDindexed_contract(1.0, opr0, shape(0,1,2), mpo0, shape(1,3,4,5), 1.0, lopr, shape(0,3,5,2,4));

    btas::QSTmergeInfo<2, Q> q_mg_bra(btas::make_array(q_l_bra, q_n_bra), btas::make_array(d_l_bra, d_n_bra));
    btas::QSTmergeInfo<2, Q> q_mg_ket(btas::make_array(q_l_ket, q_n_ket), btas::make_array(d_l_ket, d_n_ket));

    btas::QSDArray<4, Q> oprx;
    btas::QSTmerge(q_mg_bra, lopr, oprx);
    btas::QSTmerge(oprx, q_mg_ket, oprm);
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

    btas::QSDArray<5, Q> ropr;
    btas::QSDindexed_contract(1.0, mpo0, shape(0,1,2,3), opr0, shape(4,3,5), 1.0, ropr, shape(1,4,0,2,5));

    btas::QSTmergeInfo<2, Q> q_mg_bra(btas::make_array(q_n_bra, q_r_bra), btas::make_array(d_n_bra, d_r_bra));
    btas::QSTmergeInfo<2, Q> q_mg_ket(btas::make_array(q_n_ket, q_r_ket), btas::make_array(d_n_ket, d_r_ket));

    btas::QSDArray<4, Q> oprx;
    btas::QSTmerge(q_mg_bra, ropr, oprx);
    btas::QSTmerge(oprx, q_mg_ket, oprm);
  }
}

template<class Q>
double optimize_onesite_merged
(bool forward, const btas::QSDArray<4, Q>& mpo0,
                     btas::QSDArray<3, Q>& lopr,
                     btas::QSDArray<3, Q>& ropr,
                     btas::QSDArray<3, Q>& wfn0,
                     btas::QSDArray<3, Q>& wfn1,
               const double& T, const int& M = 0)
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  btas::QSDArray<3, Q> lopr_mg;
  btas::QSDArray<3, Q> ropr_mg;
  btas::QSDArray<2, Q> wfnc_mg;
  btas::QSTmergeInfo<2, Q> q_mg_ket;

  time_stamp ts;

  cout << "\t\tconstructing merged super blocks..." << flush;
  if(forward) {
    const btas::Qshapes<Q>& q_l_ket =-lopr.qshape(2);
    const btas::Qshapes<Q>& q_n_ket =-mpo0.qshape(2);
    const btas::Dshapes&    d_l_ket = lopr.dshape(0);
    const btas::Dshapes&    d_n_ket = mpo0.dshape(1);
    q_mg_ket.reset(btas::make_array(q_l_ket, q_n_ket), btas::make_array(d_l_ket, d_n_ket));

    btas::QSTmerge(q_mg_ket, wfn0, wfnc_mg);

    compute_merged_block(1, mpo0, lopr, lopr_mg);
    ropr_mg.reference(ropr);
  }
  else {
    const btas::Qshapes<Q>& q_n_ket =-mpo0.qshape(2);
    const btas::Qshapes<Q>& q_r_ket =-ropr.qshape(2);
    const btas::Dshapes&    d_n_ket = mpo0.dshape(1);
    const btas::Dshapes&    d_r_ket = ropr.dshape(0);
    q_mg_ket.reset(btas::make_array(q_n_ket, q_r_ket), btas::make_array(d_n_ket, d_r_ket));

    btas::QSTmerge(wfn0, q_mg_ket, wfnc_mg);

    lopr_mg.reference(lopr);
    compute_merged_block(0, mpo0, ropr, ropr_mg);
  }
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\tcomputing diagonal elements..." << flush;
  btas::QSDArray<2, Q> diag(wfnc_mg.q(), wfnc_mg.qshape());
  compute_diagonal_elements(lopr_mg, ropr_mg, diag);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << flush;
  davidson::Functor<2, Q> f_sigmavector = boost::bind(compute_sigmavector<Q>, lopr_mg, ropr_mg, _1, _2);
  double energy = davidson::diagonalize(f_sigmavector, diag, wfnc_mg);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    btas::QSDArray<2, Q> lmps_mg;
    btas::QSDArray<2, Q> gaug;
    canonicalize(1, wfnc_mg, lmps_mg, gaug, M);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << flush;
    btas::QSTexpand(q_mg_ket, lmps_mg, wfn0);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tcomputing guess wavefunction to the next..." << flush;
    btas::QSDArray<3, Q> wfn1_tmp;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, gaug, wfn1, 1.0, wfn1_tmp);
    wfn1 = wfn1_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    lopr.clear();
    renormalize(1, lopr_mg, lmps_mg, lmps_mg, lopr);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    btas::QSDArray<2, Q> rmps_mg;
    btas::QSDArray<2, Q> gaug;
    canonicalize(0, wfnc_mg, rmps_mg, gaug, M);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << flush;
    btas::QSTexpand(rmps_mg, q_mg_ket, wfn0);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\tcomputing guess wavefunction to the next..." << flush;
    btas::QSDArray<3, Q> wfn1_tmp;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, wfn1, gaug, 1.0, wfn1_tmp);
    wfn1 = wfn1_tmp;
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    ropr.clear();
    renormalize(0, ropr_mg, rmps_mg, rmps_mg, ropr);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

template<class Q>
double optimize_twosite_merged
(bool forward, const btas::QSDArray<4, Q>& lmpo,
               const btas::QSDArray<4, Q>& rmpo,
                     btas::QSDArray<3, Q>& lopr,
                     btas::QSDArray<3, Q>& ropr,
                     btas::QSDArray<3, Q>& lwfn,
                     btas::QSDArray<3, Q>& rwfn,
               const double& T, const int& M = 0)
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  time_stamp ts;

  cout << "\t\tconstructing merged super blocks..." << flush;

  const btas::Qshapes<Q>& q_l_ket =-lopr.qshape(2);
  const btas::Qshapes<Q>& q_m_ket =-lmpo.qshape(2);
  const btas::Dshapes&    d_l_ket = lopr.dshape(0);
  const btas::Dshapes&    d_m_ket = lmpo.dshape(1);
  btas::QSTmergeInfo<2, Q> q_lmg_ket(btas::make_array(q_l_ket, q_m_ket), btas::make_array(d_l_ket, d_m_ket));

  btas::QSDArray<3, Q> lopr_mg;
  compute_merged_block(1, lmpo, lopr, lopr_mg);

  const btas::Qshapes<Q>& q_n_ket =-rmpo.qshape(2);
  const btas::Qshapes<Q>& q_r_ket =-ropr.qshape(2);
  const btas::Dshapes&    d_n_ket = rmpo.dshape(1);
  const btas::Dshapes&    d_r_ket = ropr.dshape(0);
  btas::QSTmergeInfo<2, Q> q_rmg_ket(btas::make_array(q_n_ket, q_r_ket), btas::make_array(d_n_ket, d_r_ket));

  btas::QSDArray<3, Q> ropr_mg;
  compute_merged_block(0, rmpo, ropr, ropr_mg);

  btas::QSDArray<2, Q> wfnc_mg;
  {
    btas::QSDArray<4, Q> wfn2;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, lwfn, rwfn, 1.0, wfn2);
    btas::QSTmerge(q_lmg_ket, wfn2, q_rmg_ket, wfnc_mg);
  }

  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\tcomputing diagonal elements..." << flush;
  btas::QSDArray<2, Q> diag(wfnc_mg.q(), wfnc_mg.qshape());
  compute_diagonal_elements(lopr_mg, ropr_mg, diag);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  cout << "\t\toptimizing wavefunction (Davidson solver)..." << flush;
  davidson::Functor<2, Q> f_sigmavector = boost::bind(compute_sigmavector<Q>, lopr_mg, ropr_mg, _1, _2);
  double energy = davidson::diagonalize(f_sigmavector, diag, wfnc_mg);
  cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

  btas::QSDArray<2, Q> lwfn_mg;
  btas::QSDArray<2, Q> rwfn_mg;

  if(forward) {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    canonicalize(1, wfnc_mg, lwfn_mg, rwfn_mg, M);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << flush;
    btas::QSTexpand(q_lmg_ket, lwfn_mg, lwfn);
    btas::QSTexpand(rwfn_mg, q_rmg_ket, rwfn);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    lopr.clear();
    renormalize(1, lopr_mg, lwfn_mg, lwfn_mg, lopr);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  else {
    cout << "\t\tdoing singular value decomposition on wavefunction..." << flush;
    canonicalize(0, wfnc_mg, rwfn_mg, lwfn_mg, M);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\texpanding wavefunction..." << flush;
    btas::QSTexpand(q_lmg_ket, lwfn_mg, lwfn);
    btas::QSTexpand(rwfn_mg, q_rmg_ket, rwfn);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;

    cout << "\t\trenormalizing operators to the next..." << flush;
    ropr.clear();
    renormalize(0, ropr_mg, rwfn_mg, rwfn_mg, ropr);
    cout << "done ( " << fixed << setprecision(2) << setw(8) << ts.lap() << " sec. ) " << endl;
  }
  cout << "\t\t--------------------------------------------------------------------------------" << endl;
  cout << "\t\tTotal time for optimization: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. " << endl;

  return energy;
}

template<class Q>
double dmrg_sweep(MpOperators<Q>& mpos, MpStates<Q>& mpss, DMRGALGORITHM algo, const DmrgInput& input)
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

  assert(mpos.size() == N);
  assert(mpss.size() == N);

  btas::QSDArray<3, Q> lopr;
  btas::QSDArray<3, Q> ropr;

  double emin = 1.0e8;

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  load(mpos[0], get_mpofile(input.prefix, 0));
  load(mpss[0], get_mpsfile(input.prefix, WAVEFUNCTION,  0));
  load(lopr,    get_oprfile(input.prefix, LEFTCANONICAL, 0));

  for(size_t i = 0; i < N-1; ++i) {

    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSITE [ " << setw(3) << i << " ] " << endl;
    cout << "\t----------------------------------------------------------------------------------------------------" << endl;

    cout << "\t\tloading operators and wavefunction of next site (env)..." << flush;
    load(mpos[i+1], get_mpofile(input.prefix, i+1));
    load(mpss[i+1], get_mpsfile(input.prefix, RIGHTCANONICAL, i+1));
    cout << "done" << endl;

    // diagonalize
    double eswp;
    if(algo == ONESITE) {
      cout << "\t\toptimizing wavefunction: 1-site algorithm " << endl;
      load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, i));
//    eswp = optimize_onesite(1, mpos[i],            lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
      eswp = optimize_onesite_merged(1, mpos[i],            lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
    }
    else {
      cout << "\t\toptimizing wavefunction: 2-site algorithm " << endl;
      load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, i+1));
//    eswp = optimize_twosite(1, mpos[i], mpos[i+1], lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
      eswp = optimize_twosite_merged(1, mpos[i], mpos[i+1], lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
    }
    if(eswp < emin) emin = eswp;
    // print result
    cout << "\t\t--------------------------------------------------------------------------------" << endl;
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    cout << "\t\t--------------------------------------------------------------------------------" << endl;

    cout << "\t\tsaving operators and wavefunction of this site (sys)..." << flush;
    save(mpss[i], get_mpsfile(input.prefix, LEFTCANONICAL, i));
    save(lopr,    get_oprfile(input.prefix, LEFTCANONICAL, i+1));
    mpos[i].clear();
    mpss[i].clear();
    cout << "done" << endl;
  }
  save(mpss[N-1], get_mpsfile(input.prefix, WAVEFUNCTION, 0));

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, N-1));

  for(size_t i = N-1; i > 0; --i) {

    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSITE [ " << setw(3) << i << " ] " << endl;
    cout << "\t----------------------------------------------------------------------------------------------------" << endl;

    cout << "\t\tloading operators and wavefunction of next site (env)..." << flush;
    load(mpos[i-1], get_mpofile(input.prefix, i-1));
    load(mpss[i-1], get_mpsfile(input.prefix, LEFTCANONICAL, i-1));
    cout << "done" << endl;

    // diagonalize
    double eswp;
    if(algo == ONESITE) {
      cout << "\t\toptimizing wavefunction: 1-site algorithm " << endl;
      load(lopr,      get_oprfile(input.prefix, LEFTCANONICAL, i));
//    eswp = optimize_onesite(0, mpos[i],            lopr, ropr, mpss[i], mpss[i-1], 0.1*T, M);
      eswp = optimize_onesite_merged(0, mpos[i],            lopr, ropr, mpss[i], mpss[i-1], 0.1*T, M);
    }
    else {
      cout << "\t\toptimizing wavefunction: 2-site algorithm " << endl;
      load(lopr,      get_oprfile(input.prefix, LEFTCANONICAL, i-1));
//    eswp = optimize_twosite(0, mpos[i-1], mpos[i], lopr, ropr, mpss[i-1], mpss[i], 0.1*T, M);
      eswp = optimize_twosite_merged(0, mpos[i-1], mpos[i], lopr, ropr, mpss[i-1], mpss[i], 0.1*T, M);
    }
    if(eswp < emin) emin = eswp;
    // print result
    cout << "\t\t--------------------------------------------------------------------------------" << endl;
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    cout << "\t\t--------------------------------------------------------------------------------" << endl;

    cout << "\t\tsaving operators and wavefunction for this site (sys)..." << flush;
    save(mpss[i], get_mpsfile(input.prefix, RIGHTCANONICAL, i));
    save(ropr,    get_oprfile(input.prefix, RIGHTCANONICAL, i-1));
    mpos[i].clear();
    mpss[i].clear();
    cout << "done" << endl;
  }
  save(mpss[0], get_mpsfile(input.prefix, WAVEFUNCTION, 0));
  mpos[0].clear();
  mpss[0].clear();

  cout << "\t====================================================================================================" << endl;

  return emin;
}

double dmrg(const DmrgInput& input);

}; // namespace mpsxx

#endif // _MPSXX_CXX11_DMRG_H
