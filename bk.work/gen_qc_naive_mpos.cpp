#include <iostream>
#include <iomanip>

#include "gen_qc_naive_mpos.h"

namespace mpsxx {

qn_array_type braq;
qn_array_type ketq;

}

std::vector<int> mpsxx::gen_qc_naive_mpos (
  size_t N_sites,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint,
  std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos)
{
  braq.resize(4);
  braq[0] = fermion( 0, 0);
  braq[1] = fermion( 1, 1);
  braq[2] = fermion( 1,-1);
  braq[3] = fermion( 2, 0);

  ketq.resize(4);
  ketq[0] = fermion( 0, 0);
  ketq[1] = fermion(-1,-1);
  ketq[2] = fermion(-1, 1);
  ketq[3] = fermion(-2, 0);

  // mpos[term][site]

  size_t maxt;
  maxt  = 2*N_sites*N_sites; // 1-el terms
  maxt += N_sites;
  if(N_sites > 1) maxt +=  9*N_sites*(N_sites-1);
  if(N_sites > 2) maxt += 10*N_sites*(N_sites-1)*(N_sites-2);
  if(N_sites > 3) maxt +=  4*N_sites*(N_sites-1)*(N_sites-2)*(N_sites-3);

  mpos.clear();
  mpos.reserve(maxt);

  std::vector<int> groups;
  groups.reserve(maxt);

  gen_v1_ii(N_sites,oneint,mpos,groups);
  gen_v1_ij(N_sites,oneint,mpos,groups);

  gen_v2_iiii(N_sites,twoint,mpos,groups);
  gen_v2_iiij(N_sites,twoint,mpos,groups);
  gen_v2_ijjj(N_sites,twoint,mpos,groups);
  gen_v2_iijj(N_sites,twoint,mpos,groups);
  gen_v2_iijk(N_sites,twoint,mpos,groups);
  gen_v2_ijjk(N_sites,twoint,mpos,groups);
  gen_v2_ijkk(N_sites,twoint,mpos,groups);
  gen_v2_ijkl(N_sites,twoint,mpos,groups);

  return groups;
}

mpsxx::qn_array_type mpsxx::make_mpo_qshape (
  const mpsxx::qn_array_type& lq,
  const mpsxx::qn_array_type& iq,
        qn_shape_type& tq,
        bool last)
{
  qn_array_type iqtot;

  tq[0] = lq;
  tq[1] = braq;
  tq[2] = ketq;

  if(last) {
    tq[3] = qn_array_type(1,fermion::zero());
  }
  else {
    iqtot = lq * iq;
    tq[3] = -iqtot;
  }

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::Iden (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 0, 0));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,0,0),vp);
  mpo.insert(btas::make_array(0,1,1,0),vp);
  mpo.insert(btas::make_array(0,2,2,0),vp);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 1, 1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;
  btas::TArray<double,4> vm(1,1,1,1); vm(0,0,0,0) =-factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,0,0),vp);
  mpo.insert(btas::make_array(0,3,2,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 1,-1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,0,0),vp);
  mpo.insert(btas::make_array(0,3,1,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion(-1,-1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;
  btas::TArray<double,4> vm(1,1,1,1); vm(0,0,0,0) =-factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,1,0),vp);
  mpo.insert(btas::make_array(0,2,3,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion(-1, 1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,2,0),vp);
  mpo.insert(btas::make_array(0,1,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_CreA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 2, 0));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,0,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_CreB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  return CreB_CreA(-factor,lq,mpos,last);
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 0, 0));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,1,0),vp);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 0, 2));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,2,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 0,-2));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,1,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 0, 0));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,2,0),vp);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesA_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion(-2, 0));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesB_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  return DesA_DesB(-factor,lq,mpos,last);
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_CreB_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 1, 1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vm(1,1,1,1); vm(0,0,0,0) =-factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,2,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_CreA_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 1,-1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,1,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesA_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion(-1, 1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_DesB_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion(-1,-1));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vm(1,1,1,1); vm(0,0,0,0) =-factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,3,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesA_CreB_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<mpsxx::MPO<double,fermion>>& mpos, bool last)
{
  qn_shape_type tq;
  qn_array_type iq(1,fermion( 0, 0));
  qn_array_type iqtot = make_mpo_qshape(lq,iq,tq,last);

  dn_shape_type td;
  td[0] = dn_array_type(tq[0].size(),1);
  td[1] = dn_array_type(tq[1].size(),1);
  td[2] = dn_array_type(tq[2].size(),1);
  td[3] = dn_array_type(tq[3].size(),1);

  btas::TArray<double,4> vp(1,1,1,1); vp(0,0,0,0) = factor;

  MPO<double,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

std::vector<mpsxx::OpFunctor> mpsxx::get_op_task (mpsxx::OpFuncBase op1, size_t n, size_t i)
{
  std::vector<OpFunctor> ops;

  size_t s = 0;
  for(; s < n-1; ++s) {
    if(s == i)
      ops.push_back(boost::bind(op1,_1,_2,_3,false));
    else
      ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,false));
  }
  if(s == i)
    ops.push_back(boost::bind(op1,_1,_2,_3,true));
  else
    ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,true));

  return ops;
}

std::vector<mpsxx::OpFunctor> mpsxx::get_op_task (mpsxx::OpFuncBase op1, mpsxx::OpFuncBase op2, size_t n, size_t i, size_t j)
{
  std::vector<OpFunctor> ops;

  size_t s = 0;
  for(; s < n-1; ++s) {
    if(s == i)
      ops.push_back(boost::bind(op1,_1,_2,_3,false));
    else if(s == j)
      ops.push_back(boost::bind(op2,_1,_2,_3,false));
    else
      ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,false));
  }
  if(s == i)
    ops.push_back(boost::bind(op1,_1,_2,_3,true));
  else if(s == j)
    ops.push_back(boost::bind(op2,_1,_2,_3,true));
  else
    ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,true));

  return ops;
}

std::vector<mpsxx::OpFunctor> mpsxx::get_op_task (mpsxx::OpFuncBase op1, mpsxx::OpFuncBase op2, mpsxx::OpFuncBase op3, size_t n, size_t i, size_t j, size_t k)
{
  std::vector<OpFunctor> ops;

  size_t s = 0;
  for(; s < n-1; ++s) {
    if(s == i)
      ops.push_back(boost::bind(op1,_1,_2,_3,false));
    else if(s == j)
      ops.push_back(boost::bind(op2,_1,_2,_3,false));
    else if(s == k)
      ops.push_back(boost::bind(op3,_1,_2,_3,false));
    else
      ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,false));
  }
  if(s == i)
    ops.push_back(boost::bind(op1,_1,_2,_3,true));
  else if(s == j)
    ops.push_back(boost::bind(op2,_1,_2,_3,true));
  else if(s == k)
    ops.push_back(boost::bind(op3,_1,_2,_3,true));
  else
    ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,true));

  return ops;
}

std::vector<mpsxx::OpFunctor> mpsxx::get_op_task (mpsxx::OpFuncBase op1, mpsxx::OpFuncBase op2, mpsxx::OpFuncBase op3, mpsxx::OpFuncBase op4, size_t n, size_t i, size_t j, size_t k, size_t l)
{
  std::vector<OpFunctor> ops;

  size_t s = 0;
  for(; s < n-1; ++s) {
    if(s == i)
      ops.push_back(boost::bind(op1,_1,_2,_3,false));
    else if(s == j)
      ops.push_back(boost::bind(op2,_1,_2,_3,false));
    else if(s == k)
      ops.push_back(boost::bind(op3,_1,_2,_3,false));
    else if(s == l)
      ops.push_back(boost::bind(op4,_1,_2,_3,false));
    else
      ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,false));
  }
  if(s == i)
    ops.push_back(boost::bind(op1,_1,_2,_3,true));
  else if(s == j)
    ops.push_back(boost::bind(op2,_1,_2,_3,true));
  else if(s == k)
    ops.push_back(boost::bind(op3,_1,_2,_3,true));
  else if(s == l)
    ops.push_back(boost::bind(op4,_1,_2,_3,true));
  else
    ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,true));

  return ops;
}

int mpsxx::generate_mpos (double factor, std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos, std::vector<mpsxx::OpFunctor>& ops)
{
  if(fabs(factor) < 1.0e-16) return -1;

  size_t N_sites = ops.size();

  std::vector<MPO<double,fermion>> tmp;
  tmp.reserve(N_sites);

  qn_array_type lq(1,fermion( 0, 0));
  for(size_t i = 0; i < N_sites-1; ++i) lq = ops[i](1.0,lq,tmp);
  ops[N_sites-1](factor,lq,tmp);

  mpos.push_back(tmp);

  return 0;
}

/// group control for 1-el integrals.
/// here, all terms in 1-el part are defined as the same group
int mpsxx::get_group (size_t n, size_t i, size_t j, mpogen::OpString str)
{
  using namespace mpogen;
  size_t h = n/2;
  int g;
  if(i == j) {
    g = 0;
  }
  else {
    if((str & 0x03) == 0x03) g = 1       + j;
    if((str & 0x03) == 0x02) g = 1 +   n + j;
    if((str & 0x03) == 0x01) g = 1 + 2*n + j;
    if((str & 0x03) == 0x00) g = 1 + 3*n + j;
  }
  return g;
}

/// group control for 2-el integrals.
/// here, indices are divided by the last index 'l'
int mpsxx::get_group (size_t n, size_t i, size_t j, size_t k, size_t l, mpogen::OpString str)
{
  using namespace mpogen;
  size_t h = n/2;
  int g;
  // iiii
  if(i == j && j == k && k == l) {
    g = 0;
  }
//else {
//  if((str & 0x03) == 0x03) g = 1       + l;
//  if((str & 0x03) == 0x02) g = 1 +   n + l;
//  if((str & 0x03) == 0x01) g = 1 + 2*n + l;
//  if((str & 0x03) == 0x00) g = 1 + 3*n + l;
//}
  // iiij
  else if(i == j && j == k) {
    if((str & 0x03) == 0x03) g = 1       + l;
    if((str & 0x03) == 0x02) g = 1 +   n + l;
    if((str & 0x03) == 0x01) g = 1 + 2*n + l;
    if((str & 0x03) == 0x00) g = 1 + 3*n + l;
  }
  // ijjj
  else if(j == k && k == l) {
    if((str & 0xc0) == 0xc0) g = 1 + 4*n + i;
    if((str & 0xc0) == 0x80) g = 1 + 5*n + i;
    if((str & 0xc0) == 0x40) g = 1 + 6*n + i;
    if((str & 0xc0) == 0x00) g = 1 + 7*n + i;
  }
  // iijj
  else if(i == j && k == l) {
    if((str & 0x0f) == 0x0f) g = 1 + 8*n + k*n + l;
    if((str & 0x0f) == 0x0e) g = 1 + 8*n + k*n + l +   n*n;
    if((str & 0x0f) == 0x0d) g = 1 + 8*n + k*n + l + 2*n*n;
    if((str & 0x0f) == 0x0c) g = 1 + 8*n + k*n + l + 3*n*n;
    if((str & 0x0f) == 0x0b) g = 1 + 8*n + k*n + l + 4*n*n;
    if((str & 0x0f) == 0x0a) g = 1 + 8*n + k*n + l + 5*n*n;
    if((str & 0x0f) == 0x09) g = 1 + 8*n + k*n + l + 6*n*n;
    if((str & 0x0f) == 0x08) g = 1 + 8*n + k*n + l + 7*n*n;
    if((str & 0x0f) == 0x07) g = 1 + 8*n + k*n + l + 8*n*n;
    if((str & 0x0f) == 0x06) g = 1 + 8*n + k*n + l + 9*n*n;
    if((str & 0x0f) == 0x05) g = 1 + 8*n + k*n + l +10*n*n;
    if((str & 0x0f) == 0x04) g = 1 + 8*n + k*n + l +11*n*n;
    if((str & 0x0f) == 0x03) g = 1 + 8*n + k*n + l +12*n*n;
    if((str & 0x0f) == 0x02) g = 1 + 8*n + k*n + l +13*n*n;
    if((str & 0x0f) == 0x01) g = 1 + 8*n + k*n + l +14*n*n;
    if((str & 0x0f) == 0x00) g = 1 + 8*n + k*n + l +15*n*n;
  }
  // iijk
  else if(i == j) {
    if((str & 0x0f) == 0x0f) g = 1 + 8*n + k*n + l;
    if((str & 0x0f) == 0x0e) g = 1 + 8*n + k*n + l +   n*n;
    if((str & 0x0f) == 0x0d) g = 1 + 8*n + k*n + l + 2*n*n;
    if((str & 0x0f) == 0x0c) g = 1 + 8*n + k*n + l + 3*n*n;
    if((str & 0x0f) == 0x0b) g = 1 + 8*n + k*n + l + 4*n*n;
    if((str & 0x0f) == 0x0a) g = 1 + 8*n + k*n + l + 5*n*n;
    if((str & 0x0f) == 0x09) g = 1 + 8*n + k*n + l + 6*n*n;
    if((str & 0x0f) == 0x08) g = 1 + 8*n + k*n + l + 7*n*n;
    if((str & 0x0f) == 0x07) g = 1 + 8*n + k*n + l + 8*n*n;
    if((str & 0x0f) == 0x06) g = 1 + 8*n + k*n + l + 9*n*n;
    if((str & 0x0f) == 0x05) g = 1 + 8*n + k*n + l +10*n*n;
    if((str & 0x0f) == 0x04) g = 1 + 8*n + k*n + l +11*n*n;
    if((str & 0x0f) == 0x03) g = 1 + 8*n + k*n + l +12*n*n;
    if((str & 0x0f) == 0x02) g = 1 + 8*n + k*n + l +13*n*n;
    if((str & 0x0f) == 0x01) g = 1 + 8*n + k*n + l +14*n*n;
    if((str & 0x0f) == 0x00) g = 1 + 8*n + k*n + l +15*n*n;
  }
  // ijjk
  else if(j == k) {
    if((str & 0xc0) == 0xc0) g = 1 + 4*n + i;
    if((str & 0xc0) == 0x80) g = 1 + 5*n + i;
    if((str & 0xc0) == 0x40) g = 1 + 6*n + i;
    if((str & 0xc0) == 0x00) g = 1 + 7*n + i;
  }
  // ijkk
  else if(k == l) {
    if((str & 0x0f) == 0x0f) g = 1 + 8*n + k*n + l;
    if((str & 0x0f) == 0x0e) g = 1 + 8*n + k*n + l +   n*n;
    if((str & 0x0f) == 0x0d) g = 1 + 8*n + k*n + l + 2*n*n;
    if((str & 0x0f) == 0x0c) g = 1 + 8*n + k*n + l + 3*n*n;
    if((str & 0x0f) == 0x0b) g = 1 + 8*n + k*n + l + 4*n*n;
    if((str & 0x0f) == 0x0a) g = 1 + 8*n + k*n + l + 5*n*n;
    if((str & 0x0f) == 0x09) g = 1 + 8*n + k*n + l + 6*n*n;
    if((str & 0x0f) == 0x08) g = 1 + 8*n + k*n + l + 7*n*n;
    if((str & 0x0f) == 0x07) g = 1 + 8*n + k*n + l + 8*n*n;
    if((str & 0x0f) == 0x06) g = 1 + 8*n + k*n + l + 9*n*n;
    if((str & 0x0f) == 0x05) g = 1 + 8*n + k*n + l +10*n*n;
    if((str & 0x0f) == 0x04) g = 1 + 8*n + k*n + l +11*n*n;
    if((str & 0x0f) == 0x03) g = 1 + 8*n + k*n + l +12*n*n;
    if((str & 0x0f) == 0x02) g = 1 + 8*n + k*n + l +13*n*n;
    if((str & 0x0f) == 0x01) g = 1 + 8*n + k*n + l +14*n*n;
    if((str & 0x0f) == 0x00) g = 1 + 8*n + k*n + l +15*n*n;
  }
  // ijkl
  else {
    if((str & 0x0f) == 0x0f) g = 1 + 8*n + k*n + l;
    if((str & 0x0f) == 0x0e) g = 1 + 8*n + k*n + l +   n*n;
    if((str & 0x0f) == 0x0d) g = 1 + 8*n + k*n + l + 2*n*n;
    if((str & 0x0f) == 0x0c) g = 1 + 8*n + k*n + l + 3*n*n;
    if((str & 0x0f) == 0x0b) g = 1 + 8*n + k*n + l + 4*n*n;
    if((str & 0x0f) == 0x0a) g = 1 + 8*n + k*n + l + 5*n*n;
    if((str & 0x0f) == 0x09) g = 1 + 8*n + k*n + l + 6*n*n;
    if((str & 0x0f) == 0x08) g = 1 + 8*n + k*n + l + 7*n*n;
    if((str & 0x0f) == 0x07) g = 1 + 8*n + k*n + l + 8*n*n;
    if((str & 0x0f) == 0x06) g = 1 + 8*n + k*n + l + 9*n*n;
    if((str & 0x0f) == 0x05) g = 1 + 8*n + k*n + l +10*n*n;
    if((str & 0x0f) == 0x04) g = 1 + 8*n + k*n + l +11*n*n;
    if((str & 0x0f) == 0x03) g = 1 + 8*n + k*n + l +12*n*n;
    if((str & 0x0f) == 0x02) g = 1 + 8*n + k*n + l +13*n*n;
    if((str & 0x0f) == 0x01) g = 1 + 8*n + k*n + l +14*n*n;
    if((str & 0x0f) == 0x00) g = 1 + 8*n + k*n + l +15*n*n;
  }
  return g;
}

void mpsxx::gen_v1_ii (
  const size_t& N_sites,
  const btas::TArray<double,2>& oneint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites; ++i) {
    ops = get_op_task(mpogen::CreA_DesA,N_sites,i);
    if(generate_mpos( oneint(i,i),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,mpogen::CD));

    ops = get_op_task(mpogen::CreB_DesB,N_sites,i);
    if(generate_mpos( oneint(i,i),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,mpogen::cd));
  }
}

void mpsxx::gen_v1_ij (
  const size_t& N_sites,
  const btas::TArray<double,2>& oneint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 2) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {
      ops = get_op_task(mpogen::CreA,mpogen::DesA,N_sites,i,j);
      if(generate_mpos( oneint(i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,mpogen::CD));

      ops = get_op_task(mpogen::CreB,mpogen::DesB,N_sites,i,j);
      if(generate_mpos( oneint(i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,mpogen::cd));

      ops = get_op_task(mpogen::DesA,mpogen::CreA,N_sites,i,j);
      if(generate_mpos(-oneint(i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,mpogen::DC));

      ops = get_op_task(mpogen::DesB,mpogen::CreB,N_sites,i,j);
      if(generate_mpos(-oneint(i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,mpogen::dc));
    }
  }
}

void mpsxx::gen_v2_iiii (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites; ++i) {
    ops = get_op_task(mpogen::CreA_DesA_CreB_DesB,N_sites,i);
    if(generate_mpos( twoint(i,i,i,i),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,i,i,mpogen::CDcd));
  }
}

void mpsxx::gen_v2_iiij (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 2) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {
      ops = get_op_task(mpogen::CreA_CreB_DesB,mpogen::DesA,N_sites,i,j);
      if(generate_mpos( twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,i,j,mpogen::CcdD));

      ops = get_op_task(mpogen::CreB_CreA_DesA,mpogen::DesB,N_sites,i,j);
      if(generate_mpos( twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,i,j,mpogen::cCDd));

      ops = get_op_task(mpogen::CreA_DesA_DesB,mpogen::CreB,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,i,j,mpogen::CDdc));

      ops = get_op_task(mpogen::CreB_DesB_DesA,mpogen::CreA,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,i,j,mpogen::cdDC));
    }
  }
}

void mpsxx::gen_v2_ijjj (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 2) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {
      ops = get_op_task(mpogen::CreA,mpogen::CreB_DesB_DesA,N_sites,i,j);
      if(generate_mpos( twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,j,mpogen::CcdD));

      ops = get_op_task(mpogen::CreB,mpogen::CreA_DesA_DesB,N_sites,i,j);
      if(generate_mpos( twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,j,mpogen::cCDd));

      ops = get_op_task(mpogen::DesA,mpogen::CreA_CreB_DesB,N_sites,i,j);
      if(generate_mpos(-twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,j,mpogen::DCcd));

      ops = get_op_task(mpogen::DesB,mpogen::CreB_CreA_DesA,N_sites,i,j);
      if(generate_mpos(-twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,j,mpogen::dcCD));
    }
  }
}

void mpsxx::gen_v2_iijj (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 2) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {
      ops = get_op_task(mpogen::CreA_CreB,mpogen::DesA_DesB,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::CcDd));

      ops = get_op_task(mpogen::DesA_DesB,mpogen::CreA_CreB,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::DdCc));

      ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA_DesA,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::CDCD));

      ops = get_op_task(mpogen::CreA_DesB,mpogen::CreB_DesA,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::CdcD));

      ops = get_op_task(mpogen::CreB_DesA,mpogen::CreA_DesB,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::cDCd));

      ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB_DesB,N_sites,i,j);
      if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::cdcd));

      ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA_DesA,N_sites,i,j);
      if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::CDCD));

      ops = get_op_task(mpogen::CreA_DesA,mpogen::CreB_DesB,N_sites,i,j);
      if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::CDcd));

      ops = get_op_task(mpogen::CreB_DesB,mpogen::CreA_DesA,N_sites,i,j);
      if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::cdCD));

      ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB_DesB,N_sites,i,j);
      if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,j,mpogen::cdcd));
    }
  }
}

void mpsxx::gen_v2_iijk (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 3) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-2; ++i) {
    for(size_t j = i+1; j < N_sites-1; ++j) {
      for(size_t k = j+1; k < N_sites; ++k) {
        ops = get_op_task(mpogen::CreA_CreB,mpogen::DesA,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CcDd));

        ops = get_op_task(mpogen::CreB_CreA,mpogen::DesB,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cCdD));

        ops = get_op_task(mpogen::DesA_DesB,mpogen::CreA,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::DdCc));

        ops = get_op_task(mpogen::DesB_DesA,mpogen::CreB,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::dDcC));

        ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CDCD));

        ops = get_op_task(mpogen::CreA_DesB,mpogen::CreB,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CdcD));

        ops = get_op_task(mpogen::CreB_DesA,mpogen::CreA,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cDCd));

        ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cdcd));

        ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CDCD));

        ops = get_op_task(mpogen::CreA_DesA,mpogen::CreB,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CDcd));

        ops = get_op_task(mpogen::CreB_DesB,mpogen::CreA,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cdCD));

        ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cdcd));

        ops = get_op_task(mpogen::CreA_DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CDDC));

        ops = get_op_task(mpogen::CreA_DesB,mpogen::DesA,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CdDc));

        ops = get_op_task(mpogen::CreB_DesA,mpogen::DesB,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cDdC));

        ops = get_op_task(mpogen::CreB_DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cddc));

        ops = get_op_task(mpogen::CreA_DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CDDC));

        ops = get_op_task(mpogen::CreA_DesA,mpogen::DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::CDdc));

        ops = get_op_task(mpogen::CreB_DesB,mpogen::DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cdDC));

        ops = get_op_task(mpogen::CreB_DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,i,j,k,mpogen::cddc));
      }
    }
  }
}

void mpsxx::gen_v2_ijjk (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 3) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-2; ++i) {
    for(size_t j = i+1; j < N_sites-1; ++j) {
      for(size_t k = j+1; k < N_sites; ++k) {
        ops = get_op_task(mpogen::CreA,mpogen::DesA_DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::CDdc));

        ops = get_op_task(mpogen::CreB,mpogen::DesB_DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::cdDC));

        ops = get_op_task(mpogen::DesA,mpogen::CreA_CreB,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::DCcd));

        ops = get_op_task(mpogen::DesB,mpogen::CreB_CreA,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::dcCD));

        ops = get_op_task(mpogen::CreA,mpogen::CreA_DesA,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::CCDD));

        ops = get_op_task(mpogen::CreA,mpogen::CreB_DesA,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::CcDd));

        ops = get_op_task(mpogen::CreB,mpogen::CreA_DesB,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::cCdD));

        ops = get_op_task(mpogen::CreB,mpogen::CreB_DesB,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::ccdd));

        ops = get_op_task(mpogen::CreA,mpogen::CreA_DesA,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::CCDD));

        ops = get_op_task(mpogen::CreA,mpogen::CreB_DesB,mpogen::DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::CcdD));

        ops = get_op_task(mpogen::CreB,mpogen::CreA_DesA,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::cCDd));

        ops = get_op_task(mpogen::CreB,mpogen::CreB_DesB,mpogen::DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::ccdd));

        ops = get_op_task(mpogen::DesA,mpogen::CreA_DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::DCDC));

        ops = get_op_task(mpogen::DesA,mpogen::CreA_DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::DCdc));

        ops = get_op_task(mpogen::DesB,mpogen::CreB_DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::dcDC));

        ops = get_op_task(mpogen::DesB,mpogen::CreB_DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::dcdc));

        ops = get_op_task(mpogen::DesA,mpogen::CreA_DesA,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::DCDC));

        ops = get_op_task(mpogen::DesA,mpogen::CreB_DesB,mpogen::CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::DcdC));

        ops = get_op_task(mpogen::DesB,mpogen::CreA_DesA,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::dCDc));

        ops = get_op_task(mpogen::DesB,mpogen::CreB_DesB,mpogen::CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,j,k,mpogen::dcdc));
      }
    }
  }
}

void mpsxx::gen_v2_ijkk (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 3) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-2; ++i) {
    for(size_t j = i+1; j < N_sites-1; ++j) {
      for(size_t k = j+1; k < N_sites; ++k) {
        ops = get_op_task(mpogen::CreA,mpogen::CreB,mpogen::DesA_DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::CcDd));

        ops = get_op_task(mpogen::CreB,mpogen::CreA,mpogen::DesB_DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::cCdD));

        ops = get_op_task(mpogen::DesA,mpogen::DesB,mpogen::CreA_CreB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::DdCc));

        ops = get_op_task(mpogen::DesB,mpogen::DesA,mpogen::CreB_CreA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::dDcC));

        ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA_DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::CDCD));

        ops = get_op_task(mpogen::CreA,mpogen::DesB,mpogen::CreB_DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::CdcD));

        ops = get_op_task(mpogen::CreB,mpogen::DesA,mpogen::CreA_DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::cDCd));

        ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB_DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::cdcd));

        ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA_DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::CDCD));

        ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreB_DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::CDcd));

        ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreA_DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::cdCD));

        ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB_DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::cdcd));

        ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA_DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::DCCD));

        ops = get_op_task(mpogen::DesA,mpogen::CreB,mpogen::CreA_DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::DcCd));

        ops = get_op_task(mpogen::DesB,mpogen::CreA,mpogen::CreB_DesA,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::dCcD));

        ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB_DesB,N_sites,i,j,k);
        if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::dccd));

        ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA_DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::DCCD));

        ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreB_DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::DCcd));

        ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreA_DesA,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::dcCD));

        ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB_DesB,N_sites,i,j,k);
        if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,k,mpogen::dccd));
      }
    }
  }
}

void mpsxx::gen_v2_ijkl (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<int>& groups)
{
  if(N_sites < 4) return;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-3; ++i) {
    for(size_t j = i+1; j < N_sites-2; ++j) {
      for(size_t k = j+1; k < N_sites-1; ++k) {
        for(size_t l = k+1; l < N_sites; ++l) {
          ops = get_op_task(mpogen::CreA,mpogen::CreA,mpogen::DesA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CCDD));

          ops = get_op_task(mpogen::CreA,mpogen::CreB,mpogen::DesA,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CcDd));

          ops = get_op_task(mpogen::CreB,mpogen::CreA,mpogen::DesB,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cCdD));

          ops = get_op_task(mpogen::CreB,mpogen::CreB,mpogen::DesB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::ccdd));

          ops = get_op_task(mpogen::CreA,mpogen::CreA,mpogen::DesA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CCDD));

          ops = get_op_task(mpogen::CreA,mpogen::CreB,mpogen::DesB,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CcdD));

          ops = get_op_task(mpogen::CreB,mpogen::CreA,mpogen::DesA,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cCDd));

          ops = get_op_task(mpogen::CreB,mpogen::CreB,mpogen::DesB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::ccdd));

          ops = get_op_task(mpogen::DesA,mpogen::DesA,mpogen::CreA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DDCC));

          ops = get_op_task(mpogen::DesA,mpogen::DesB,mpogen::CreA,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DdCc));

          ops = get_op_task(mpogen::DesB,mpogen::DesA,mpogen::CreB,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dDcC));

          ops = get_op_task(mpogen::DesB,mpogen::DesB,mpogen::CreB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::ddcc));

          ops = get_op_task(mpogen::DesA,mpogen::DesA,mpogen::CreA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DDCC));

          ops = get_op_task(mpogen::DesA,mpogen::DesB,mpogen::CreB,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DdcC));

          ops = get_op_task(mpogen::DesB,mpogen::DesA,mpogen::CreA,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dDCc));

          ops = get_op_task(mpogen::DesB,mpogen::DesB,mpogen::CreB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::ddcc));

          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CDDC));

          ops = get_op_task(mpogen::CreA,mpogen::DesB,mpogen::DesA,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CdDc));

          ops = get_op_task(mpogen::CreB,mpogen::DesA,mpogen::DesB,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cDdC));

          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cddc));

          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CDDC));

          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CDdc));

          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cdDC));

          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cddc));

          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DCCD));

          ops = get_op_task(mpogen::DesA,mpogen::CreB,mpogen::CreA,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DcCd));

          ops = get_op_task(mpogen::DesB,mpogen::CreA,mpogen::CreB,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dCcD));

          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dccd));

          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DCCD));

          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DCcd));

          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dcCD));

          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dccd));

          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CDCD));

          ops = get_op_task(mpogen::CreA,mpogen::DesB,mpogen::CreB,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CdcD));

          ops = get_op_task(mpogen::CreB,mpogen::DesA,mpogen::CreA,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cDCd));

          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cdcd));

          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CDCD));

          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::CDcd));

          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cdCD));

          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::cdcd));

          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DCDC));

          ops = get_op_task(mpogen::DesA,mpogen::CreB,mpogen::DesB,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DcdC));

          ops = get_op_task(mpogen::DesB,mpogen::CreA,mpogen::DesA,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dCDc));

          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dcdc));

          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DCDC));

          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::DCdc));

          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dcDC));

          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
          if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(get_group(N_sites,i,j,k,l,mpogen::dcdc));
        }
      }
    }
  }
}
