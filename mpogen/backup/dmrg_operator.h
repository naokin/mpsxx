#ifndef _MPSXX_DMRG_OPERATOR_H
#define _MPSXX_DMRG_OPERATOR_H 1

namespace mpsxx {

enum OpType
{
  ZERO        ,
  IDENTITY    ,
  HAM         ,
  CRE         ,
  DES         ,
  CRE_COMP    ,
  DES_COMP    ,
  CRE_CRE     ,
  CRE_DES     ,
  DES_CRE     ,
  DES_DES     ,
  CRE_CRE_COMP,
  CRE_DES_COMP,
  DES_CRE_COMP,
  DES_DES_COMP 
};

template<OpType L_optype = ZERO>
struct OpComponent {
  template<OpType R_Op> constexpr static OpType site_optype() { return ZERO; }
};

template<>
struct OpComponent<IDENTITY> {
  template<OpType R_Op> constexpr static OpType site_optype() { return R_Op; }
};

template<>
struct OpComponent<HAM> {
  template<> constexpr static OpType site_optype<IDENTITY>() { return IDENTITY; }
};

template<>
struct OpComponent<CRE> {
  template<> constexpr static OpType site_optype<CRE>() { return IDENTITY; }
  template<> constexpr static OpType site_optype<HAM>() { return DES_COMP; }
};

template<>
struct OpComponent<DES> {
  template<> constexpr static OpType site_optype<DES>() { return IDENTITY; }
  template<> constexpr static OpType site_optype<HAM>() { return CRE_COMP; }
};

template<OpType L_Op, OpType R_Op>
constexpr OpType get_site_optype() { return OpComponent<L_Op>::site_optype<R_Op>(); }

}; // namespace mpsxx

#endif // _MPSXX_DMRG_OPERATOR_H
