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
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, wfn0, wfn1, 1.0, wfn2);

    diag.resize(wfn2.q(), wfn2.qshape());
    compute_diagonal_elements(mpo0, mpo1, lopr, ropr, diag);

    f_sigmavector = boost::bind(compute_sigmavector<Q>, mpo0, mpo1, lopr, ropr, _1, _2);
  }
  else {
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, wfn1, wfn0, 1.0, wfn2);

    diag.resize(wfn2.q(), wfn2.qshape());
    compute_diagonal_elements(mpo1, mpo0, lopr, ropr, diag);

    f_sigmavector = boost::bind(compute_sigmavector<Q>, mpo1, mpo0, lopr, ropr, _1, _2);
  }

  double energy = davidson::diagonalize(f_sigmavector, diag, wfn2);

  if(forward) {
    canonicalize(1, wfn2, wfn0, wfn1, M);

    btas::QSDArray<3, Q> lopr_tmp;
    renormalize(1, mpo0, lopr, wfn0, wfn0, lopr_tmp);

    lopr = lopr_tmp;
  }
  else {
    canonicalize(0, wfn2, wfn0, wfn1, M);

    btas::QSDArray<3, Q> ropr_tmp;
    renormalize(0, mpo0, ropr, wfn0, wfn0, ropr_tmp);

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

  load(mpos[0], get_mpofile(input.prefix, MOLECULAR,     0));
  load(mpss[0], get_mpsfile(input.prefix, WAVEFUNCTION,  0));
  load(lopr,    get_oprfile(input.prefix, LEFTCANONICAL, 0));

  for(size_t i = 0; i < N-1; ++i) {

    load(mpos[i+1], get_mpofile(input.prefix, MOLECULAR,      i));
    load(mpss[i+1], get_mpsfile(input.prefix, RIGHTCANONICAL, i));
    load(ropr,      get_oprfile(input.prefix, RIGHTCANONICAL, i));

    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(1, mpos[i],            lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
    else                eswp = optimize_twosite(1, mpos[i], mpos[i+1], lopr, ropr, mpss[i], mpss[i+1], 0.1*T, M);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;

    save(mpss[i], get_mpsfile(input.prefix, LEFTCANONICAL, i));
    save(lopr,    get_oprfile(input.prefix, LEFTCANONICAL, i+1));
    mpos[i].clear();
    mpss[i].clear();
  }
  save(mpss[N-1], get_mpsfile(input.prefix, WAVEFUNCTION, 0));

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  load(ropr, get_oprfile(input.prefix, RIGHTCANONICAL, N-1));

  for(size_t i = N-1; i > 0; --i) {

    load(mpos[i-1], get_mpofile(input.prefix, MOLECULAR,     i));
    load(mpss[i-1], get_mpsfile(input.prefix, LEFTCANONICAL, i));
    load(lopr,      get_oprfile(input.prefix, LEFTCANONICAL, i));

    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(0, mpos[i],            lopr, ropr, mpss[i], mpss[i-1], 0.1*T, M);
    else                eswp = optimize_twosite(0, mpos[i], mpos[i-1], lopr, ropr, mpss[i], mpss[i-1], 0.1*T, M);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;

    save(mpss[i], get_mpsfile(input.prefix, RIGHTCANONICAL, i));
    save(ropr,    get_oprfile(input.prefix, RIGHTCANONICAL, i-1));
    mpos[i].clear();
    mpss[i].clear();
  }
  save(mpss[0], get_mpsfile(input.prefix, WAVEFUNCTION, 0));
  mpos[0].clear();
  mpss[0].clear();

  return emin;
}

double dmrg(const DmrgInput& input);

}; // namespace mpsxx

#endif // _MPSXX_CXX11_DMRG_H
