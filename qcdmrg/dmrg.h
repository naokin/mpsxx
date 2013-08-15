#ifndef _MPSXX_CXX11_DMRG_H
#define _MPSXX_CXX11_DMRG_H 1

#include <iostream>
#include <iomanip>

#include <btas/QSPARSE/QSDArray.h>

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
  btas::QSDArray<3, Q> diag(wfn0.q(), wfn0.qshape());
  compute_diagonal_elements(mpo0, lopr, ropr, diag);

  davidson::Functor<3, Q> f_sigmavector = boost::bind(compute_sigmavector<Q>, mpo0, lopr, ropr, _1, _2);
  double energy = davidson::diagonalize(f_sigmavector, diag, wfn0);

  if(forward) {
    btas::QSDArray<3, Q> lmps;
    canonicalize(1, wfn0, lmps, M);

    btas::QSDArray<3, Q> wfn1_tmp;
    compute_guesswave(1, lmps, wfn0, wfn1, wfn1_tmp);

    btas::QSDArray<3, Q> lopr_tmp;
    renormalize(1, mpo0, lopr, lmps, lmps, lopr_tmp);

    wfn0 = lmps;
    wfn1 = wfn1_tmp;
    lopr = lopr_tmp;
  }
  else {
    btas::QSDArray<3, Q> rmps;
    canonicalize(0, wfn0, rmps, M);

    btas::QSDArray<3, Q> wfn1_tmp;
    compute_guesswave(0, rmps, wfn0, wfn1, wfn1_tmp);

    btas::QSDArray<3, Q> ropr_tmp;
    renormalize(0, mpo0, ropr, rmps, rmps, ropr_tmp);

    wfn0 = rmps;
    wfn1 = wfn1_tmp;
    ropr = ropr_tmp;
  }

  return energy;
}

template<class Q>
double optimize_twosite
(bool forward, const btas::QSDArray<4, Q>& mpo0,
               const btas::QSDArray<4, Q>& mpo1,
                     btas::QSDArray<3, Q>& lopr,
                     btas::QSDArray<3, Q>& ropr,
                     btas::QSDArray<3, Q>& wfn0,
                     btas::QSDArray<3, Q>& wfn1,
               const double& T, const int& M = 0)
{
  btas::QSDArray<4, Q> wfn2;
  btas::QSDArray<4, Q> diag;

  davidson::Functor<4, Q> f_sigmavector;

  if(forward) {
//MPSXX_DEBUG("DEBUG[optimize_twosite]::01");
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, wfn0, wfn1, 1.0, wfn2);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::02");
    diag.resize(wfn2.q(), wfn2.qshape());
    compute_diagonal_elements(mpo0, mpo1, lopr, ropr, diag);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::03");
    f_sigmavector = boost::bind(compute_sigmavector<Q>, mpo0, mpo1, lopr, ropr, _1, _2);
  }
  else {
//MPSXX_DEBUG("DEBUG[optimize_twosite]::04");
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, wfn1, wfn0, 1.0, wfn2);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::05");
    diag.resize(wfn2.q(), wfn2.qshape());
    compute_diagonal_elements(mpo1, mpo0, lopr, ropr, diag);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::06");
    f_sigmavector = boost::bind(compute_sigmavector<Q>, mpo1, mpo0, lopr, ropr, _1, _2);
  }

//MPSXX_DEBUG("DEBUG[optimize_twosite]::07");
  double energy = davidson::diagonalize(f_sigmavector, diag, wfn2);

  if(forward) {
//MPSXX_DEBUG("DEBUG[optimize_twosite]::08");
    canonicalize(1, wfn2, wfn0, wfn1, M);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::09");
    btas::QSDArray<3, Q> lopr_tmp;
    renormalize(1, mpo0, lopr, wfn0, wfn0, lopr_tmp);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::10");
    lopr = lopr_tmp;
  }
  else {
//MPSXX_DEBUG("DEBUG[optimize_twosite]::11");
    canonicalize(0, wfn2, wfn0, wfn1, M);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::12");
    btas::QSDArray<3, Q> ropr_tmp;
    renormalize(0, mpo0, ropr, wfn0, wfn0, ropr_tmp);

//MPSXX_DEBUG("DEBUG[optimize_twosite]::13");
    ropr = ropr_tmp;
  }

  return energy;
}

template<class Q>
double dmrg_sweep(MpOperators<Q>& mpos, MpStates<Q>& mpss, DMRGALGORITHM algo, const DmrgInput& input)
{
  using std::cout;
  using std::endl;
  using std::setw;
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

//MPSXX_DEBUG("DEBUG[dmrg_sweep]::01");
//load(mpos[0], get_mpofile(input.prefix, MOLECULAR,     0));
  load(mpos[0], get_mpofile(input.prefix, 0));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::02");
  load(mpss[0], get_mpsfile(input.prefix, WAVEFUNCTION,  0));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::03");
  load(lopr,    get_oprfile(input.prefix, LEFTCANONICAL, 0));

  for(size_t i = 0; i < N-1; ++i) {

//MPSXX_DEBUG("DEBUG[dmrg_sweep]::04");
//  load(mpos[i+1], get_mpofile(input.prefix, MOLECULAR,      i+1));
    load(mpos[i+1], get_mpofile(input.prefix, i+1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::05");
    load(mpss[i+1], get_mpsfile(input.prefix, RIGHTCANONICAL, i+1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::06");

    // diagonalize
    double eswp;
    if(algo == ONESITE) {
      load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, i));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::07");
      eswp = optimize_onesite(1, mpos[i],            lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
    }
    else {
      load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, i+1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::07");
      eswp = optimize_twosite(1, mpos[i], mpos[i+1], lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
    }
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::08");
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;

//MPSXX_DEBUG("DEBUG[dmrg_sweep]::09");
    save(mpss[i], get_mpsfile(input.prefix, LEFTCANONICAL, i));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::10");
    save(lopr,    get_oprfile(input.prefix, LEFTCANONICAL, i+1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::11");
    mpos[i].clear();
    mpss[i].clear();
  }
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::12");
  save(mpss[N-1], get_mpsfile(input.prefix, WAVEFUNCTION, 0));

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

//MPSXX_DEBUG("DEBUG[dmrg_sweep]::13");
  load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, N-1));

  for(size_t i = N-1; i > 0; --i) {

//MPSXX_DEBUG("DEBUG[dmrg_sweep]::14");
//  load(mpos[i-1], get_mpofile(input.prefix, MOLECULAR,     i-1));
    load(mpos[i-1], get_mpofile(input.prefix, i-1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::15");
    load(mpss[i-1], get_mpsfile(input.prefix, LEFTCANONICAL, i-1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::16");

    // diagonalize
    double eswp;
    if(algo == ONESITE) {
      load(lopr,      get_oprfile(input.prefix, LEFTCANONICAL, i));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::17");
      eswp = optimize_onesite(0, mpos[i],            lopr, ropr, mpss[i], mpss[i-1], 0.1*T, M);
    }
    else {
      load(lopr,      get_oprfile(input.prefix, LEFTCANONICAL, i-1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::17");
      eswp = optimize_twosite(0, mpos[i], mpos[i-1], lopr, ropr, mpss[i], mpss[i-1], 0.1*T, M);
    }
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::18");
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;

//MPSXX_DEBUG("DEBUG[dmrg_sweep]::19");
    save(mpss[i], get_mpsfile(input.prefix, RIGHTCANONICAL, i));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::20");
    save(ropr,    get_oprfile(input.prefix, RIGHTCANONICAL, i-1));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::21");
    mpos[i].clear();
    mpss[i].clear();
  }
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::22");
  save(mpss[0], get_mpsfile(input.prefix, WAVEFUNCTION, 0));
//MPSXX_DEBUG("DEBUG[dmrg_sweep]::23");
  mpos[0].clear();
  mpss[0].clear();

  return emin;
}

double dmrg(const DmrgInput& input);

}; // namespace mpsxx

#endif // _MPSXX_CXX11_DMRG_H
