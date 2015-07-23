#ifndef __MPSXX_GEN_QC_NAIVE_MPOS_H
#define __MPSXX_GEN_QC_NAIVE_MPOS_H

#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <MPX_Types.h>
#include <symmetry/fermion.h>

namespace mpsxx {

typedef btas::Qshapes<fermion> qn_array_type;

typedef btas::TVector<qn_array_type,4> qn_shape_type;

typedef btas::Dshapes dn_array_type;

typedef btas::TVector<dn_array_type,4> dn_shape_type;

typedef boost::function<qn_array_type(const double&,const qn_array_type&,std::vector<MPO<double,fermion>>&,bool)> OpFuncBase;

typedef boost::function<qn_array_type(const double&,const qn_array_type&,std::vector<MPO<double,fermion>>&)> OpFunctor;

std::vector<int> gen_qc_naive_mpos (size_t N_sites, const btas::TArray<double,2>& oneint, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos);

qn_array_type make_mpo_qshape (const qn_array_type& lq, const qn_array_type& iq, qn_shape_type& tq, bool last);

namespace mpogen {

enum OpString {
  CD=0xd,
  cd=0x8,
  DC=0x7,
  dc=0x2,

  CCDD=0xf5,
  CcDd=0xe4,
  cCdD=0xb1,
  CcdD=0xe1,
  cCDd=0xb4,
  ccdd=0xa0,

  CDCD=0xdd,
  CDcd=0xd8,
  cdCD=0x8d,
  CdcD=0xc9,
  cDCd=0x9c,
  cdcd=0x88,

  CDDC=0xd7,
  CDdc=0xd2,
  cdDC=0x87,
  CdDc=0xc6,
  cDdC=0x93,
  cddc=0x82,

  DCCD=0x7d,
  DCcd=0x78,
  dcCD=0x2d,
  DcCd=0x6c,
  dCcD=0x39,
  dccd=0x28,

  DCDC=0x77,
  DCdc=0x72,
  dcDC=0x27,
  DcdC=0x63,
  dCDc=0x36,
  dcdc=0x22,

  DDCC=0x5f,
  DdCc=0x4e,
  dDcC=0x1b,
  DdcC=0x4b,
  dDCc=0x1e,
  ddcc=0x0a
};

qn_array_type Iden (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type DesA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);

qn_array_type CreB_CreA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreA_CreB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreA_DesA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreA_DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreB_DesA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreB_DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type DesA_DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type DesB_DesA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);

qn_array_type CreA_CreB_DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreB_CreA_DesA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreA_DesA_DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);
qn_array_type CreB_DesB_DesA (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);

qn_array_type CreA_DesA_CreB_DesB (const double& factor, const qn_array_type& lq, std::vector<MPO<double,fermion>>& mpos, bool last);

} // namespace mpogen

std::vector<OpFunctor> get_op_task (OpFuncBase op1, size_t n, size_t i);
std::vector<OpFunctor> get_op_task (OpFuncBase op1, OpFuncBase op2, size_t n, size_t i, size_t j);
std::vector<OpFunctor> get_op_task (OpFuncBase op1, OpFuncBase op2, OpFuncBase op3, size_t n, size_t i, size_t j, size_t k);
std::vector<OpFunctor> get_op_task (OpFuncBase op1, OpFuncBase op2, OpFuncBase op3, OpFuncBase op4, size_t n, size_t i, size_t j, size_t k, size_t l);

int generate_mpos (double factor, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<OpFunctor>& ops);

/// group control for 1-el integrals.
/// here, all terms in 1-el part are defined as the same group
int get_group (size_t n, size_t i, size_t j, mpogen::OpString str);

/// group control for 2-el integrals.
/// here, indices are divided by the last index 'l'
int get_group (size_t n, size_t i, size_t j, size_t k, size_t l, mpogen::OpString str);

void gen_v1_ii (const size_t& N_sites, const btas::TArray<double,2>& oneint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v1_ij (const size_t& N_sites, const btas::TArray<double,2>& oneint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);

void gen_v2_iiii (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_iiij (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_ijjj (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_iijj (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_iijk (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_ijjk (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_ijkk (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);
void gen_v2_ijkl (const size_t& N_sites, const btas::TArray<double,4>& twoint, std::vector<std::vector<MPO<double,fermion>>>& mpos, std::vector<int>& groups);

} // namespace mpsxx

#endif // __MPSXX_GEN_QC_NAIVE_MPOS_H
