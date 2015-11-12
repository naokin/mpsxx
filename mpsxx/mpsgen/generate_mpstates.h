#ifndef __MPSXX_GENERATE_MPSTATES_H
#define __MPSXX_GENERATE_MPSTATES_H 1

#include <MpOperators.h>
#include <MpStates.h>
#include <MpSite.h>

namespace mpsxx {

/// Remove non-contributed quantum numbers for each boundary
template<class Q>
void remove_zero_boundary (
      const std::vector<btas::Qshapes<Q>>& qn, ///< physical indices
            std::vector<btas::Qshapes<Q>>& qb, ///< boundary indices
            size_t __Max_quanta = 0)
{
   // reduce zero quantum blocks with backward direction
   for(size_t i = L-1; i > 0; --i)
   {
      btas::Qshapes<Q>& qx = qn[i];
      btas::Qshapes<Q>& ql = qb[i-1];
      btas::Qshapes<Q>& qr = qb[i];

      // check non-zero for each ql index
      auto lt = ql.begin();
      while(lt != ql.end())
      {
         bool non_zero = false;
         for(size_t p = 0; p < qx.size(); ++p)
            for(size_t r = 0; r < qr.size(); ++r)
               non_zero |= (qr[r] == (qx[p] * (*lt)));

         if(non_zero)
            ++lt;
         else
            ql.erase(lt);
      }
      assert(ql.size() > 0);
   }

   // reduce zero quantum blocks with forward direction
   for(size_t i = 0; i < L-1; ++i)
   {
      btas::Qshapes<Q>& qx = qn[i];
      btas::Qshapes<Q>& ql = qb[i];
      btas::Qshapes<Q>& qr = qb[i+1];

      // further reduction
      if(__Max_quanta > 0 && ql.size() > __Max_quanta)
      {
         size_t offs = (ql.size()-__Max_quanta)/2;
         ql = btas::Qshapes<Q>(ql.begin()+offs, ql.begin()+offs+__Max_quanta);

         // check non-zero for each qr index
         auto rt = qr.begin();
         while(rt != qr.end())
         {
            bool non_zero = false;
            for(size_t l = 0; l < ql.size(); ++l)
               for(size_t p = 0; p < qx.size(); ++p)
                  non_zero |= (*rt == (ql[l] * qx[p]));

            if(non_zero)
               ++rt;
            else
               qr.erase(rt);
         }
         assert(qr.size() > 0);
      }
   }
}

/// Generate quantum numbers for each boundary
/// \param L lattice size
/// \param qt total quantum number of lattice
/// \param __Max_quanta maximum number of quantum blocks on each boundary (0: no limitation by default)
template<class Q>
std::vector<btas::Qshapes<Q>> generate_quanta (const size_t& L, const Q& qt, size_t __Max_quanta = 0)
{
   // physical indices
   std::vector<btas::Qshapes<Q>> qn(L, MpSite<Q>::quanta());

   // boundary indices
   std::vector<btas::Qshapes<Q>> qb(L);

   // generate quantum number blocks for MPS
   qb[0] = qn[0];
   for(size_t i = 1; i < L-1; ++i)
   {
      // get unique elements of { q(l) x q(n) }
      qb[i] = qb[i-1] & qn[i];
   }
   qb[L-1] = btas::Qshapes<Q>(1, qt); // used for checking quantum states

   remove_zero_boundary(qn, qb, __Max_quanta);

   return qb;
}

/// Generate quantum numbers for each boundary from MPO
/// \param mpos MPO's on memory
/// \param qt total quantum number of lattice
/// \param __Max_quanta maximum number of quantum blocks on each boundary (0: no limitation by default)
template<class Q>
std::vector<btas::Qshapes<Q>> generate_quantum_states (const MpOperator<Q>& mpos, const Q& qt, size_t __Max_quanta = 0)
{
   size_t L = mpos.size();

   // physical indices
   std::vector<btas::Qshapes<Q>> qn(L);

   for(size_t i = 0; i < L; ++i)
   {
      qn[i] = -mpo[i].qshape(2);
   }

   // boundary quantum numbers
   std::vector<btas::Qshapes<Q>> qb(L);

   // generate quantum number blocks for MPS
   qb[0] = qn[0];
   for(size_t i = 1; i < L-1; ++i)
   {
      // get unique elements of { q(l) x q(n) }
      qb[i] = qb[i-1] & qn[i];
   }
   qb[L-1] = btas::Qshapes<Q>(1, qt); // used for checking quantum states

   remove_zero_boundary(qn, qb, __Max_quanta);

   return qb;
}

/// Generate MPS by Generator
/// \param prefix file prefix to save MPS objects
/// \param L lattice size
/// \param qt total quantum number of lattice
/// \param gen e.g. random number generator called by double gen();
/// \param __Max_quanta maximum number of quantum blocks on each boundary (0: no limitation by default)
template<class Q, class Generator>
void generate_mpstates_on_disk (const std::string& prefix, const size_t& L, const Q& qt, Generator gen, size_t __Max_quanta = 0)
{
   btas::Qshapes<Q> qn(MpSite<Q>::quanta());
   btas::Dshapes    dn(qn.size(), 1);

   btas::Qshapes<Q> ql(1, Q::zero());
   btas::Dshapes    dl(ql.size(), 1);

   btas::Qshapes<Q> qr;
   btas::Dshapes    dr;

   std::vector<btas::Qshapes<Q>> qb = generate_quanta(L, qt, __Max_quanta);

   btas::TVector<btas::Qshapes<Q>, 3> q_shape;
   btas::TVector<btas::Dshapes,    3> d_shape;

   for(size_t i = 0; i < L-1; ++i)
   {
      qr = qb[i];
      dr = btas::Dshapes(qr.size(), 1);

      q_shape = btas::make_array( ql, qn,-qr);
      d_shape = btas::make_array( dl, dn,-dr);

      btas::QSDArray<3, Q> mps(Q::zero(), q_shape, d_shape, gen);
      save(mps, get_mpsfile(prefix, LEFTCANONICAL, i));

      ql = qr;
      dl = dr;
   }

   // Last site
   {
      qr = btas::Qshapes<Q>(1, Q::zero()); // qb[L-1] is not used
      dr = btas::Dshapes   (qr.size(), 1);

      q_shape = btas::make_array( ql, qn,-qr);
      d_shape = btas::make_array( dl, dn,-dr);

      btas::QSDArray<3, Q> mps(qt, q_shape, d_shape, gen);
      btas::Normalize(mps);
      save(mps, get_mpsfile(prefix, WAVEFUNCTION, L-1));
   }

   btas::QSDArray<3, Q> lmps;
   btas::QSDArray<3, Q> rmps;

   load(rmps, get_mpsfile(prefix, WAVEFUNCTION, L-1));

   // initial canonicalization
   for(size_t i = L-1; i > 0; --i)
   {
      btas::QSDArray<2, Q> gauge;
      canonicalize(0, rmps, gauge, __Max_quanta);
      save(rmps, get_mpsfile(prefix, RIGHTCANONICAL, i));

      load(lmps, get_mpsfile(prefix, LEFTCANONICAL, i-1));

      btas::QSDArray<3, Q> wfnc;
      btas::Contract(1.0, lmps, btas::shape(2), gauge, btas::shape(0), 1.0, wfnc);
      btas::Normalize(wfnc);

      rmps = wfnc;
   }

   save(rmps, get_mpsfile(prefix, WAVEFUNCTION, 0));
}

//! Generate MPS by Generator
/*!
 *  \param mpss MPS container
 *  \param L lattice size
 *  \param qt total quantum number of lattice
 *  \param gen e.g. random number generator called by double gen();
 *  \param __Max_quanta maximum number of quantum blocks on each boundary (0: no limitation by default)
 */
template<class Q, class Generator>
void generate_mpstates
(MpStates<Q>& mpss, const MpOperator<Q>& mpos, const Q& qt, Generator gen, size_t __Max_quanta = 0)
{
  size_t L = mpos.size();

  mpss.resize(L);

  btas::Qshapes<Q> qn(MpSite<Q>::quanta());
  btas::Dshapes    dn(qn.size(), 1);

  btas::Qshapes<Q> ql(1, Q::zero());
  btas::Dshapes    dl(ql.size(), 1);

  btas::Qshapes<Q> qr;
  btas::Dshapes    dr;

  std::vector<btas::Qshapes<Q>> qb = generate_quantum_states(L, qt, __Max_quanta);

  btas::TVector<btas::Qshapes<Q>, 3> q_shape;
  btas::TVector<btas::Dshapes,    3> d_shape;
  for(size_t i = 0; i < L-1; ++i) {
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
  qr = btas::Qshapes<Q>(1, Q::zero()); // qb[L-1] is not used
  dr = btas::Dshapes   (qr.size(), 1);

  q_shape = btas::make_array( ql, qn,-qr);
  d_shape = btas::make_array( dl, dn,-dr);

  mpss[L-1].resize(qt, q_shape, d_shape); // set qt as a total quantum number of array
  mpss[L-1].generate(gen);
  btas::QSDnormalize(mpss[L-1]);

  // initial canonicalization
  for(size_t i = L-1; i > 0; --i) {
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
