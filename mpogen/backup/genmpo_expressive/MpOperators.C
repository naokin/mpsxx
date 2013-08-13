//
// Collective Classes and Functions for Matrix Product Operators (MPO)
//

#include <MpOperators.h>



/*

template<class Q = Quantum>
btas::QSDArray<4, Q> mpsxx::MpGenerator<Q>::XXX(const OpInfo& l_op, const OpInfo& s_op, const OpInfo& r_op, ...);

*/

template<>
btas::QSDArray<4, mpsxx::Fermion> mpsxx::MpGenerator<mpsxx::Fermion>::Cre
(const OpInfo& l_op, const OpInfo& s_op, const OpInfo& r_op, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint)
{
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_cre_a
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(+1,+1));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block = 1.0; op.insert(0, site_type::alpha,  site_type::vacuum, 0, block);
  block =-1.0; op.insert(0, site_type::pair,   site_type::beta,   0, block);
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_cre_b
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(+1,-1));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block = 1.0; op.insert(0, site_type::beta,   site_type::vacuum, 0, block);
  block = 1.0; op.insert(0, site_type::pair,   site_type::alpha,  0, block);
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_des_a
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(-1,-1));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block = 1.0; op.insert(0, site_type::vacuum, site_type::alpha,  0, block);
  block =-1.0; op.insert(0, site_type::beta,   site_type::pair,   0, block);
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_des_b
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(-1,+1));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block = 1.0; op.insert(0, site_type::vacuum, site_type::beta,   0, block);
  block = 1.0; op.insert(0, site_type::alpha,  site_type::pair,   0, block);
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_cre_a_cre_b
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(+2, 0));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block =-1.0; op.insert(0, site_type::pair,   site_type::vacuum, 0, block);
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_cre_b_cre_a
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(+2, 0));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block = 1.0; op.insert(0, site_type::pair,   site_type::vacuum, 0, block);
}

btas::QSDArray<4, mpsxx::Fermion> mpsxx::_fermion_cre_a_des_a
(const mpsxx::Fermion& l_quanta, const mpsxx::Fermion& r_quanta)
{
  typedef mpsxx::MpSite<mpsxx::Fermion> site_type;
  assert((r_quanta - l_quanta) != Fermion(+2, 0));

  btas::Qshapes<Fermion> lq = { l_quanta };
  btas::Qshapes<Fermion> rq = { r_quanta };
  btas::Qshapes<Fermion> sq = Fermion::fock();

  btas::QSDArray<4, Fermion> op(Fermion::zero(), btas::make_array(lq, sq,-sq,-rq));
  btas::DArray<4> block(1,1,1,1);
  block = 1.0; op.insert(0, site_type::pair,   site_type::vacuum, 0, block);
}

