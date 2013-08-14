#ifndef _MPSXX_CXX11_GENERATE_MPSTATES_H
#define _MPSXX_CXX11_GENERATE_MPSTATES_H 1

#include <MpOperators.h>
#include <MpStates.h>
#include <MpSite.h>

namespace mpsxx {

//! Generate quantum numbers for each boundary
/*!
 *  \param N lattice size
 *  \param qt total quantum number of lattice
 *  \param _max_quantum_blocks maximum number of quantum blocks on each boundary (0: no limitation by default)
 */
template<class Q>
std::vector<btas::Qshapes<Q>> generate_quantum_states
(const size_t& N, const Q& qt, size_t _max_quantum_blocks = 0)
{
  // physical index
  btas::Qshapes<Q> qn = MpSite<Q>::quanta();
  // zero quantum number
  btas::Qshapes<Q> qz(1, Q::zero());

  // boundary quantum numbers
  std::vector<btas::Qshapes<Q>> qb(N);

  // generate quantum number blocks for MPS
  qb[0] = qn;
  for(size_t i = 1; i < N-1; ++i) {
    qb[i] = qb[i-1] & qn; // get unique elements of { q(l) x q(n) }
  }
  qb[N-1] = btas::Qshapes<Q>(1, qt); // used for checking quantum states

  // reduce zero quantum blocks
  for(size_t i = N-1; i > 0; --i) {
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
void generate_mpstates
(MpStates<Q>& mpss, const size_t& N, const Q& qt, Generator gen, size_t _max_quantum_blocks = 0)
{
  mpss.resize(N);

  btas::Qshapes<Q> qn(MpSite<Q>::quanta());
  btas::Dshapes    dn(qn.size(), 1);

  btas::Qshapes<Q> ql(1, Q::zero());
  btas::Dshapes    dl(ql.size(), 1);

  btas::Qshapes<Q> qr;
  btas::Dshapes    dr;

  std::vector<btas::Qshapes<Q>> qb = generate_quantum_states(N, qt, _max_quantum_blocks);

  btas::TVector<btas::Qshapes<Q>, 3> q_shape;
  btas::TVector<btas::Dshapes,    3> d_shape;
  for(size_t i = 0; i < N-1; ++i) {
    qr = qb[i];
    dr = btas::Dshapes(qr.size(), 1);

    q_shape = btas::make_array( ql, qn,-qr);
    d_shape = btas::make_array( dl, dn,-dr);

    mpss[i].resize(Q::zero(), q_shape, d_shape);
    mpss[i].generate(gen);
//  save(mpss[i], get_mpsfile(prefix, i));
//  mpss[i].clear();

    ql = qr;
    dl = dr;
  }
  qr = btas::Qshapes<Q>(1, Q::zero()); // qb[N-1] is not used
  dr = btas::Dshapes   (qr.size(), 1);

  q_shape = btas::make_array( ql, qn,-qr);
  d_shape = btas::make_array( dl, dn,-dr);

  mpss[N-1].resize(qt, q_shape, d_shape); // set qt as a total quantum number of array
  mpss[N-1].generate(gen);
  btas::QSDnormalize(mpss[N-1]);

  // initial canonicalization
  for(size_t i = N-1; i > 0; --i) {
    btas::QSDArray<2> gauge;
    canonicalize(0, mpss[i], gauge, M);
//  save(mpss[i], get_mpsfile(prefix, i));
//  mpss[i].clear();

//  load(mpss[i-1], get_mpsfile(prefix, i-1));
    btas::QSDArray<3> wfunc;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, mpss[i-1], gauge, 1.0, wfunc);
    btas::QSDnormalize(wfunc);
    mpss[i-1] = wfunc;
  }
//save(mpss[0], get_mpsfile(prefix, 0));
//mpss[0].clear();
}



//! Generate quantum numbers for each boundary from MPO
/*!
 *  \param N lattice size
 *  \param qt total quantum number of lattice
 *  \param _max_quantum_blocks maximum number of quantum blocks on each boundary (0: no limitation by default)
 */
template<class Q>
std::vector<btas::Qshapes<Q>> generate_quantum_states
(const MpOperator<Q>& mpos, const Q& qt, size_t _max_quantum_blocks = 0)
{
  size_t N = mpos.size();

  // zero quantum number
  btas::Qshapes<Q> qz(1, Q::zero());

  // boundary quantum numbers
  std::vector<btas::Qshapes<Q>> qb(N);

  // generate quantum number blocks for MPS
  qb[0] =-mpos[0].qshape(2);
  for(size_t i = 1; i < N-1; ++i) {
    qb[i] = qb[i-1] & (-mpos[i].qshape(2)); // get unique elements of { q(l) x q(n) }
  }
  qb[N-1] = btas::Qshapes<Q>(1, qt); // used for checking quantum states

  // reduce zero quantum blocks
  for(size_t i = N-1; i > 0; --i) {
    btas::Qshapes<Q>& qn =-mpos[i].qshape(2);
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
    btas::Qshapes<Q>& qn =-mpos[i].qshape(2);
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
void generate_mpstates
(MpStates<Q>& mpss, const MpOperator<Q>& mpos, const Q& qt, Generator gen, size_t _max_quantum_blocks = 0)
{
  size_t N = mpos.size();

  mpss.resize(N);

  btas::Qshapes<Q> qn(MpSite<Q>::quanta());
  btas::Dshapes    dn(qn.size(), 1);

  btas::Qshapes<Q> ql(1, Q::zero());
  btas::Dshapes    dl(ql.size(), 1);

  btas::Qshapes<Q> qr;
  btas::Dshapes    dr;

  std::vector<btas::Qshapes<Q>> qb = generate_quantum_states(N, qt, _max_quantum_blocks);

  btas::TVector<btas::Qshapes<Q>, 3> q_shape;
  btas::TVector<btas::Dshapes,    3> d_shape;
  for(size_t i = 0; i < N-1; ++i) {
    qr = qb[i];
    dr = btas::Dshapes(qr.size(), 1);

    q_shape = btas::make_array( ql, qn,-qr);
    d_shape = btas::make_array( dl, dn,-dr);

    mpss[i].resize(Q::zero(), q_shape, d_shape);
    mpss[i].generate(gen);
    save(mpss[i], get_mpsfile(prefix, i));

    mpss[i].clear();

    ql = qr;
    dl = dr;
  }
  qr = btas::Qshapes<Q>(1, Q::zero()); // qb[N-1] is not used
  dr = btas::Dshapes   (qr.size(), 1);

  q_shape = btas::make_array( ql, qn,-qr);
  d_shape = btas::make_array( dl, dn,-dr);

  mpss[N-1].resize(qt, q_shape, d_shape); // set qt as a total quantum number of array
  mpss[N-1].generate(gen);
  btas::QSDnormalize(mpss[N-1]);

  // initial canonicalization
  for(size_t i = N-1; i > 0; --i) {
    btas::QSDArray<2> gauge;
    canonicalize(0, mpss[i], gauge, M);
    save(mpss[i], get_mpsfile(prefix, i));
    mpss[i].clear();

    load(mpss[i-1], get_mpsfile(prefix, i-1));
    btas::QSDArray<3> wfunc;
    btas::QSDgemm(btas::NoTrans, btas::NoTrans, 1.0, mpss[i-1], gauge, 1.0, wfunc);
    btas::QSDnormalize(wfunc);
    mpss[i-1] = wfunc;
  }
  save(mpss[0], get_mpsfile(prefix, 0));
  mpss[0].clear();
}

};

#endif // _MPSXX_CXX11_GENERATE_MPSTATES_H
