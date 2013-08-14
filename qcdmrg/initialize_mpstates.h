#ifndef _MPSXX_CXX11_INITIALIZE_MPSTATES_H
#define _MPSXX_CXX11_INITIALIZE_MPSTATES_H 1

#include <MpOperators.h>
#include <MpStates.h>

#include <driver/canonicalize.h>
#include <driver/renormalize.h>
#include <driver/guesswave.h>
#include <driver/fileio.h>

namespace mpsxx {

//! Generate quantum numbers for each boundary from MPO
/*!
 *  \param sites quantum numbers for each site (physical index)
 *  \param qt total quantum number of lattice
 *  \param _max_quantum_blocks maximum number of quantum blocks on each boundary (0: no limitation by default)
 */
template<class Q>
std::vector<btas::Qshapes<Q>> generate_quantum_states
(const std::vector<btas::Qshapes<Q>>& sites, const Q& qt, size_t _max_quantum_blocks = 0)
{
  size_t N = sites.size();

  // zero quantum number
  btas::Qshapes<Q> qz(1, Q::zero());

  // boundary quantum numbers
  std::vector<btas::Qshapes<Q>> qb(N);

  // generate quantum number blocks for MPS
  qb[0] = sites[0];
  for(size_t i = 1; i < N-1; ++i) {
    qb[i] = qb[i-1] & sites[i]; // get unique elements of { q(l) x q(n) }
  }
  qb[N-1] = btas::Qshapes<Q>(1, qt); // used for checking quantum states

  // reduce zero quantum blocks
  for(size_t i = N-1; i > 0; --i) {
    const btas::Qshapes<Q>& qn = sites[i];
          btas::Qshapes<Q>& ql = qb[i-1];
          btas::Qshapes<Q>& qr = qb[i];

    // check non-zero for each ql index
    auto lt = ql.begin();
    while(lt != ql.end()) {
      bool non_zero = false;
      for(size_t p = 0; p < qn.size(); ++p) {
        for(size_t r = 0; r < qr.size(); ++r) {
          non_zero |= (qr[r] == (qn[p] * (*lt)));
        }
      }
      if(non_zero)
        ++lt;
      else
        ql.erase(lt);
    }
    assert(ql.size() > 0);
  }
  for(size_t i = 0; i < N-1; ++i) {
    const btas::Qshapes<Q>& qn = sites[i];
          btas::Qshapes<Q>& ql = qb[i];
          btas::Qshapes<Q>& qr = qb[i+1];

    // further reduction
    if(_max_quantum_blocks > 0 && ql.size() > _max_quantum_blocks) {
      size_t offs = (ql.size()-_max_quantum_blocks)/2;
      ql = btas::Qshapes<Q>(ql.begin()+offs, ql.begin()+offs+_max_quantum_blocks);

      // check non-zero for each qr index
      auto rt = qr.begin();
      while(rt != qr.end()) {
        bool non_zero = false;
        for(size_t l = 0; l < ql.size(); ++l) {
          for(size_t p = 0; p < qn.size(); ++p) {
            non_zero |= (*rt == (ql[l] * qn[p]));
          }
        }
        if(non_zero)
          ++rt;
        else
          qr.erase(rt);
      }
      assert(qr.size() > 0);
    }
  }
  return std::move(qb);
}

//! Generate MPS by Generator
/*!
 *  \param mpss MPS container
 *  \param N lattice size
 *  \param qt total quantum number of lattice
 *  \param gen e.g. random number generator called by double gen();
 *  \param _max_quantum_blocks maximum number of quantum blocks on each boundary (0: no limitation by default)
 */
template<class Q, class Generator>
void initialize_mpstates
(MpOperators<Q>& mpos, MpStates<Q>& mpss, const Q& qt, Generator gen, const std::string& prefix = "./", size_t _max_quantum_blocks = 0)
{
  size_t N = mpos.size();

  MPSXX_DEBUG("DEBUG[initialize_mpstates]::01");
  mpss.resize(N);

  btas::Qshapes<Q> qz(1, Q::zero());

  btas::Qshapes<Q> ql(1, Q::zero());
  btas::Dshapes    dl(ql.size(), 1);

  btas::Qshapes<Q> qr;
  btas::Dshapes    dr;

  MPSXX_DEBUG("DEBUG[initialize_mpstates]::02");
  std::vector<btas::Qshapes<Q>> qn(N);
  std::vector<btas::Dshapes   > dn(N);
  for(size_t i = 0; i < N; ++i) {
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::03");
    load(mpos[i], get_mpofile(prefix, MOLECULAR, i));
    qn[i] = -mpos[i].qshape(2);
    dn[i] = btas::Dshapes(qn[i].size(), 1);
    mpos[i].clear();
  }
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::04");
  std::vector<btas::Qshapes<Q>> qb = generate_quantum_states(qn, qt, _max_quantum_blocks);

  btas::TVector<btas::Qshapes<Q>, 3> q_shape;
  btas::TVector<btas::Dshapes,    3> d_shape;
  for(size_t i = 0; i < N-1; ++i) {
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::05");
    qr = qb[i];
    dr = btas::Dshapes(qr.size(), 1);

    q_shape = btas::make_array( ql, qn[i],-qr);
    d_shape = btas::make_array( dl, dn[i], dr);

  MPSXX_DEBUG("DEBUG[initialize_mpstates]::06");
    mpss[i].resize(Q::zero(), q_shape, d_shape);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::07");
    mpss[i].generate(gen);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::08");
    save(mpss[i], get_mpsfile(prefix, WAVEFUNCTION, i));

    mpss[i].clear();

    ql = qr;
    dl = dr;
  }
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::09");
  qr = qz; // qb[N-1] is not used
  dr = btas::Dshapes(qr.size(), 1);

  q_shape = btas::make_array( ql, qn[N-1],-qr);
  d_shape = btas::make_array( dl, dn[N-1], dr);

  MPSXX_DEBUG("DEBUG[initialize_mpstates]::10");
  mpss[N-1].resize(qt, q_shape, d_shape); // set qt as a total quantum number of array
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::11");
  mpss[N-1].generate(gen);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::12");
  btas::QSDnormalize(mpss[N-1]);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::13");
  save(mpss[N-1], get_mpsfile(prefix, WAVEFUNCTION, N-1));

  // initial canonicalization
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::14");
  btas::QSDArray<3, Q> ropr_0(Q::zero(), btas::make_array(qz, qz, qz));
  ropr_0.insert(btas::shape(0, 0, 0), btas::DArray<3>(1, 1, 1));
  ropr_0.fill(1.0);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::15");
  save(ropr_0, get_oprfile(prefix, RIGHTCANONICAL, N-1));
  for(size_t i = N-1; i > 0; --i) {
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::16");
    btas::QSDArray<3, Q> rmps;
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::17");
    canonicalize(0, mpss[i], rmps, _max_quantum_blocks);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::18");
    save(rmps, get_mpsfile(prefix, RIGHTCANONICAL, i));

    btas::QSDArray<3, Q> ropr_1;
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::19");
    load(mpos[i], get_mpofile(prefix, MOLECULAR, i));
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::20");
    renormalize(0, mpos[i], ropr_0, rmps, rmps, ropr_1);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::21");
    save(ropr_1, get_oprfile(prefix, RIGHTCANONICAL, i-1));
    ropr_0 = ropr_1;

    mpos[i].clear();

  MPSXX_DEBUG("DEBUG[initialize_mpstates]::22");
    load(mpss[i-1], get_mpsfile(prefix, WAVEFUNCTION, i-1));
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::23");
    btas::QSDArray<3, Q> lmps(mpss[i-1]);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::24");
    mpss[i-1].clear();
    compute_guesswave(0, rmps, mpss[i], lmps, mpss[i-1]);
  MPSXX_DEBUG("DEBUG[initialize_mpstates]::25");
    save(mpss[i-1], get_mpsfile(prefix, WAVEFUNCTION, i-1));

    mpss[i].clear();
  }
  mpss[0].clear();

  MPSXX_DEBUG("DEBUG[initialize_mpstates]::26");
  btas::QSDArray<3, Q> lopr_0(Q::zero(), btas::make_array(qz, qz, qz));
  lopr_0.insert(btas::shape(0, 0, 0), btas::DArray<3>(1, 1, 1));
  lopr_0.fill(1.0);
  save(lopr_0, get_oprfile(prefix, LEFTCANONICAL, 0));
}

};

#endif // _MPSXX_CXX11_INITIALIZE_MPSTATES_H
