#include <iostream>
#include <iomanip>

#include "mpidefs.h"
#include "gen_qc_naive_mpos.h"

namespace mpsxx {

qn_array_type braq;
qn_array_type ketq;

}

std::vector<int> mpsxx::gen_qc_naive_mpos (
  const size_t& N_sites,
  const double& E_core,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint,
  std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos)
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

  maxt /= Communicator().size();

  mpos.clear();
  mpos.reserve(maxt);

  std::vector<int> groups;
  groups.reserve(maxt);

  std::map<int,int> gmap;

  gen_ecore(N_sites,E_core,mpos,groups,gmap);

  gen_v1_ii(N_sites,oneint,mpos,groups,gmap);
  gen_v1_ij(N_sites,oneint,mpos,groups,gmap);

  gen_v2_iiii(N_sites,twoint,mpos,groups,gmap);
  gen_v2_iiij(N_sites,twoint,mpos,groups,gmap);
  gen_v2_ijjj(N_sites,twoint,mpos,groups,gmap);
  gen_v2_iijj(N_sites,twoint,mpos,groups,gmap);
  gen_v2_iijk(N_sites,twoint,mpos,groups,gmap);
  gen_v2_ijjk(N_sites,twoint,mpos,groups,gmap);
  gen_v2_ijkk(N_sites,twoint,mpos,groups,gmap);
  gen_v2_ijkl(N_sites,twoint,mpos,groups,gmap);

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

mpsxx::qn_array_type mpsxx::mpogen::Iden (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
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

mpsxx::qn_array_type mpsxx::mpogen::CreA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,0,0),vp);
  mpo.insert(btas::make_array(0,3,2,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,0,0),vp);
  mpo.insert(btas::make_array(0,3,1,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,1,0),vp);
  mpo.insert(btas::make_array(0,2,3,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,2,0),vp);
  mpo.insert(btas::make_array(0,1,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_CreA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,0,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_CreB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
{
  return CreB_CreA(-factor,lq,mpos,last);
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,1,0),vp);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,2,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,1,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,2,0),vp);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesA_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,0,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::DesB_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
{
  return DesA_DesB(-factor,lq,mpos,last);
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_CreB_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,2,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_CreA_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,1,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesA_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,1,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreB_DesB_DesA (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,2,3,0),vm);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

mpsxx::qn_array_type mpsxx::mpogen::CreA_DesA_CreB_DesB (const double& factor, const mpsxx::qn_array_type& lq, std::vector<btas::QSTArray<double,4,fermion>>& mpos, bool last)
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

  btas::QSTArray<double,4,fermion> mpo(fermion::zero(),tq,td,false);
  mpo.insert(btas::make_array(0,3,3,0),vp);

  std::vector<int> idx1(1,0);
  std::vector<int> idx2(1,1);
  mpo.parity(idx1,idx2);

  mpos.push_back(mpo);

  return iqtot;
}

std::vector<mpsxx::OpFunctor> mpsxx::get_op_task (size_t n)
{
  std::vector<OpFunctor> ops;

  for(size_t s = 0; s < n-1; ++s) {
    ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,false));
  }
  ops.push_back(boost::bind(mpogen::Iden,_1,_2,_3,true));

  return ops;
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

int mpsxx::generate_mpos (double factor, std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos, std::vector<mpsxx::OpFunctor>& ops)
{
  if(fabs(factor) < 1.0e-16) return -1;

  size_t N_sites = ops.size();

  std::vector<btas::QSTArray<double,4,fermion>> tmp;
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
  int m = n-1;
  int g;
  if(i == j) {
    g = 0;
  }
  else {
    if((str & 0x03) == 0x03) g =     + j;
    if((str & 0x03) == 0x02) g =   m + j;
    if((str & 0x03) == 0x01) g = 2*m + j;
    if((str & 0x03) == 0x00) g = 3*m + j;
  }
  return g;
}

/// group control for 2-el integrals.
/// here, indices are divided by the last index 'l'
int mpsxx::get_group (size_t n, size_t i, size_t j, size_t k, size_t l, mpogen::OpString str)
{
  using namespace mpogen;
  int m = n-1;
  int g;
  if(1) {
  // k-MPOs
    // iiii
    if(i == j && j == k && k == l) {
      g = 0;
    }
    // ijkl
    else {
      if((str & 0x03) == 0x03) g =     + l;
      if((str & 0x03) == 0x02) g =   m + l;
      if((str & 0x03) == 0x01) g = 2*m + l;
      if((str & 0x03) == 0x00) g = 3*m + l;
    }
  }
  else {
  // k2-MPOs
    // iiii
    if(i == j && j == k && k == l) {
      g = 0;
    }
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
  }
  return g;
}

void mpsxx::gen_ecore (
  const size_t& N_sites,
  const double& E_core,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  if(E_core == 0.0) return;

  int gid = 0;
  auto it = gmap.find(gid);

  if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
  if(it->second == iproc) {
    std::vector<OpFunctor> ops = get_op_task(N_sites);
    if(generate_mpos(E_core,mpos,ops) == 0) groups.push_back(gid);
  }
}

void mpsxx::gen_v1_ii (
  const size_t& N_sites,
  const btas::TArray<double,2>& oneint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites; ++i) {

    if(fabs(oneint(i,i)) < 1.0e-16) continue;

    gid = get_group(N_sites,i,i,mpogen::CD);
    it = gmap.find(gid);
    if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
    if(it->second == iproc) {
      ops = get_op_task(mpogen::CreA_DesA,N_sites,i);
      if(generate_mpos( oneint(i,i),mpos,ops) == 0) groups.push_back(gid);
    }

    gid = get_group(N_sites,i,i,mpogen::cd);
    it = gmap.find(gid);
    if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
    if(it->second == iproc) {
      ops = get_op_task(mpogen::CreB_DesB,N_sites,i);
      if(generate_mpos( oneint(i,i),mpos,ops) == 0) groups.push_back(gid);
    }
  }
}

void mpsxx::gen_v1_ij (
  const size_t& N_sites,
  const btas::TArray<double,2>& oneint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 2) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {

      if(fabs(oneint(i,j)) < 1.0e-16) continue;

      gid = get_group(N_sites,i,j,mpogen::CD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA,mpogen::DesA,N_sites,i,j);
        if(generate_mpos( oneint(i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,j,mpogen::cd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB,mpogen::DesB,N_sites,i,j);
        if(generate_mpos( oneint(i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,j,mpogen::DC);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::DesA,mpogen::CreA,N_sites,i,j);
        if(generate_mpos(-oneint(i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,j,mpogen::dc);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::DesB,mpogen::CreB,N_sites,i,j);
        if(generate_mpos(-oneint(i,j),mpos,ops) == 0) groups.push_back(gid);
      }
    }
  }
}

void mpsxx::gen_v2_iiii (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites; ++i) {

    if(fabs(twoint(i,i,i,i)) < 1.0e-16) continue;

    gid = get_group(N_sites,i,i,i,i,mpogen::CDcd);
    it = gmap.find(gid);
    if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
    if(it->second == iproc) {
      ops = get_op_task(mpogen::CreA_DesA_CreB_DesB,N_sites,i);
      if(generate_mpos( twoint(i,i,i,i),mpos,ops) == 0) groups.push_back(gid);
    }
  }
}

void mpsxx::gen_v2_iiij (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 2) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {

      if(fabs(twoint(i,i,i,j)) < 1.0e-16) continue;

      gid = get_group(N_sites,i,i,i,j,mpogen::CcdD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_CreB_DesB,mpogen::DesA,N_sites,i,j);
        if(generate_mpos( twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,i,j,mpogen::cCDd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB_CreA_DesA,mpogen::DesB,N_sites,i,j);
        if(generate_mpos( twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,i,j,mpogen::CDdc);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_DesA_DesB,mpogen::CreB,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,i,j,mpogen::cdDC);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB_DesB_DesA,mpogen::CreA,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,i,j),mpos,ops) == 0) groups.push_back(gid);
      }
    }
  }
}

void mpsxx::gen_v2_ijjj (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 2) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {

      if(fabs(twoint(i,j,j,j)) < 1.0e-16) continue;

      gid = get_group(N_sites,i,j,j,j,mpogen::CcdD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA,mpogen::CreB_DesB_DesA,N_sites,i,j);
        if(generate_mpos( twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,j,j,j,mpogen::cCDd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB,mpogen::CreA_DesA_DesB,N_sites,i,j);
        if(generate_mpos( twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,j,j,j,mpogen::DCcd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::DesA,mpogen::CreA_CreB_DesB,N_sites,i,j);
        if(generate_mpos(-twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,j,j,j,mpogen::dcCD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::DesB,mpogen::CreB_CreA_DesA,N_sites,i,j);
        if(generate_mpos(-twoint(i,j,j,j),mpos,ops) == 0) groups.push_back(gid);
      }
    }
  }
}

void mpsxx::gen_v2_iijj (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 2) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-1; ++i) {
    for(size_t j = i+1; j < N_sites; ++j) {

      if(fabs(twoint(i,i,j,j)) >= 1.0e-16) {

      gid = get_group(N_sites,i,i,j,j,mpogen::CcDd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_CreB,mpogen::DesA_DesB,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::DdCc);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::DesA_DesB,mpogen::CreA_CreB,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::CDCD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA_DesA,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::CdcD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_DesB,mpogen::CreB_DesA,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::cDCd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB_DesA,mpogen::CreA_DesB,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::cdcd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB_DesB,N_sites,i,j);
        if(generate_mpos(-twoint(i,i,j,j),mpos,ops) == 0) groups.push_back(gid);
      }

      } if(fabs(twoint(i,j,i,j)) >= 1.0e-16) {

      gid = get_group(N_sites,i,i,j,j,mpogen::CDCD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA_DesA,N_sites,i,j);
        if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::CDcd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreA_DesA,mpogen::CreB_DesB,N_sites,i,j);
        if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::cdCD);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB_DesB,mpogen::CreA_DesA,N_sites,i,j);
        if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      gid = get_group(N_sites,i,i,j,j,mpogen::cdcd);
      it = gmap.find(gid);
      if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
      if(it->second == iproc) {
        ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB_DesB,N_sites,i,j);
        if(generate_mpos( twoint(i,j,i,j),mpos,ops) == 0) groups.push_back(gid);
      }

      }
    }
  }
}

void mpsxx::gen_v2_iijk (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 3) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-2; ++i) {
    for(size_t j = i+1; j < N_sites-1; ++j) {
      for(size_t k = j+1; k < N_sites; ++k) {

        if(fabs(twoint(i,i,j,k)) >= 1.0e-16) {

        gid = get_group(N_sites,i,i,j,k,mpogen::CcDd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_CreB,mpogen::DesA,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cCdD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_CreA,mpogen::DesB,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::DdCc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA_DesB,mpogen::CreA,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::dDcC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB_DesA,mpogen::CreB,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CDCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CdcD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesB,mpogen::CreB,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cDCd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesA,mpogen::CreA,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cdcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CDDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CdDc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesB,mpogen::DesA,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cDdC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesA,mpogen::DesB,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cddc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,i,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        } if(fabs(twoint(i,j,i,k)) >= 1.0e-16) {

        gid = get_group(N_sites,i,i,j,k,mpogen::CDCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CDcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesA,mpogen::CreB,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cdCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesB,mpogen::CreA,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cdcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CDDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::CDdc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA_DesA,mpogen::DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cdDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesB,mpogen::DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,i,j,k,mpogen::cddc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB_DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,i,k),mpos,ops) == 0) groups.push_back(gid);
        }

        }
      }
    }
  }
}

void mpsxx::gen_v2_ijjk (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 3) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-2; ++i) {
    for(size_t j = i+1; j < N_sites-1; ++j) {
      for(size_t k = j+1; k < N_sites; ++k) {

        if(fabs(twoint(i,j,j,k)) >= 1.0e-16) {

        gid = get_group(N_sites,i,j,j,k,mpogen::CDdc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::DesA_DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::cdDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::DesB_DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::DCcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA_CreB,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::dcCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB_CreA,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::CCDD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::CreA_DesA,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::CcDd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::CreB_DesA,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::cCdD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::CreA_DesB,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::ccdd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::CreB_DesB,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::DCDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA_DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::DCdc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA_DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::dcDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB_DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::dcdc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB_DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        } if(fabs(twoint(i,j,k,j)) >= 1.0e-16) {

        gid = get_group(N_sites,i,j,j,k,mpogen::CCDD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::CreA_DesA,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::CcdD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::CreB_DesB,mpogen::DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::cCDd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::CreA_DesA,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::ccdd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::CreB_DesB,mpogen::DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::DCDC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA_DesA,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::DcdC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreB_DesB,mpogen::CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::dCDc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreA_DesA,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,j,k,mpogen::dcdc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB_DesB,mpogen::CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,j),mpos,ops) == 0) groups.push_back(gid);
        }

        }
      }
    }
  }
}

void mpsxx::gen_v2_ijkk (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 3) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-2; ++i) {
    for(size_t j = i+1; j < N_sites-1; ++j) {
      for(size_t k = j+1; k < N_sites; ++k) {

        if(fabs(twoint(i,j,k,k)) >= 1.0e-16) {

        gid = get_group(N_sites,i,j,k,k,mpogen::CcDd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::CreB,mpogen::DesA_DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::cCdD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::CreA,mpogen::DesB_DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::DdCc);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::DesB,mpogen::CreA_CreB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::dDcC);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::DesA,mpogen::CreB_CreA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::CDCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA_DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::CdcD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::DesB,mpogen::CreB_DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::cDCd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::DesA,mpogen::CreA_DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::cdcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB_DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::DCCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA_DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::DcCd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreB,mpogen::CreA_DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::dCcD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreA,mpogen::CreB_DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::dccd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB_DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,j,k,k),mpos,ops) == 0) groups.push_back(gid);
        }

        } if(fabs(twoint(i,k,j,k)) >= 1.0e-16) {

        gid = get_group(N_sites,i,j,k,k,mpogen::CDCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA_DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::CDcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreB_DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::cdCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreA_DesA,N_sites,i,j,k);
          if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::cdcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB_DesB,N_sites,i,j,k);
          if(generate_mpos( twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::DCCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA_DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::DCcd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreB_DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::dcCD);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreA_DesA,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        gid = get_group(N_sites,i,j,k,k,mpogen::dccd);
        it = gmap.find(gid);
        if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
        if(it->second == iproc) {
          ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB_DesB,N_sites,i,j,k);
          if(generate_mpos(-twoint(i,k,j,k),mpos,ops) == 0) groups.push_back(gid);
        }

        }
      }
    }
  }
}

void mpsxx::gen_v2_ijkl (
  const size_t& N_sites,
  const btas::TArray<double,4>& twoint,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<int>& groups, std::map<int,int>& gmap)
{
  if(N_sites < 4) return;

  Communicator world;
  size_t nproc = world.size();
  size_t iproc = world.rank();

  int gid;
  std::map<int,int>::iterator it;

  std::vector<OpFunctor> ops;
  for(size_t i = 0; i < N_sites-3; ++i) {
    for(size_t j = i+1; j < N_sites-2; ++j) {
      for(size_t k = j+1; k < N_sites-1; ++k) {
        for(size_t l = k+1; l < N_sites; ++l) {

          if(fabs(twoint(i,j,k,l)) >= 1.0e-16) {

          gid = get_group(N_sites,i,j,k,l,mpogen::CCDD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::CreA,mpogen::DesA,mpogen::DesA,N_sites,i,j,k,l);
              if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CcDd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::CreB,mpogen::DesA,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cCdD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::CreA,mpogen::DesB,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::ccdd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::CreB,mpogen::DesB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DDCC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::DesA,mpogen::CreA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DdCc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::DesB,mpogen::CreA,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dDcC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::DesA,mpogen::CreB,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::ddcc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::DesB,mpogen::CreB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CDDC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CdDc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesB,mpogen::DesA,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cDdC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesA,mpogen::DesB,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cddc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DCCD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DcCd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreB,mpogen::CreA,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dCcD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreA,mpogen::CreB,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dccd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,k,l),mpos,ops) == 0) groups.push_back(gid);
          }

          } if(fabs(twoint(i,j,l,k)) >= 1.0e-16) {

          gid = get_group(N_sites,i,j,k,l,mpogen::CCDD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::CreA,mpogen::DesA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CcdD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::CreB,mpogen::DesB,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cCDd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::CreA,mpogen::DesA,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::ccdd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::CreB,mpogen::DesB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DDCC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::DesA,mpogen::CreA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DdcC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::DesB,mpogen::CreB,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dDCc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::DesA,mpogen::CreA,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::ddcc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::DesB,mpogen::CreB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CDCD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CdcD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesB,mpogen::CreB,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cDCd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesA,mpogen::CreA,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cdcd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DCDC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DcdC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreB,mpogen::DesB,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dCDc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreA,mpogen::DesA,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dcdc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,j,l,k),mpos,ops) == 0) groups.push_back(gid);
          }

          } if(fabs(twoint(i,k,j,l)) >= 1.0e-16) {

          gid = get_group(N_sites,i,j,k,l,mpogen::CDDC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CDdc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cdDC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cddc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DCCD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DCcd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dcCD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dccd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos(-twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CDCD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::CDcd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreA,mpogen::DesA,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cdCD);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreA,mpogen::DesA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::cdcd);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::CreB,mpogen::DesB,mpogen::CreB,mpogen::DesB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DCDC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::DCdc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesA,mpogen::CreA,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dcDC);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::DesA,mpogen::CreA,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          gid = get_group(N_sites,i,j,k,l,mpogen::dcdc);
          it = gmap.find(gid);
          if(it == gmap.end()) it = gmap.insert(it,std::make_pair(gid,gid%nproc));
          if(it->second == iproc) {
            ops = get_op_task(mpogen::DesB,mpogen::CreB,mpogen::DesB,mpogen::CreB,N_sites,i,j,k,l);
            if(generate_mpos( twoint(i,k,j,l),mpos,ops) == 0) groups.push_back(gid);
          }

          }
        }
      }
    }
  }
}
