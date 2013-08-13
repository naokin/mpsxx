#include <iostream>
#include <cstring>

enum OpType {
  NORMAL   = 0x7fffffff,
  MASK     = 0xff000000,
  COMP     = 0x80000000,
  ZERO     = 0x00000000,
  IDN      = 0x01000000,
  HAM      = 0x02000000,
  CRE      = 0x03000000,
  DES      = 0x04000000,
  CRE_COMP = 0x83000000,
  DES_COMP = 0x84000000,
};

template<OpType T>
struct OpInfo {
  const static OpType Type = T;
  template<OpType R> struct SiteOp { const static OpType Type = ZERO; };
};

template<> template<> struct OpInfo<IDN>::SiteOp<IDN>           { const static OpType Type = IDN; };
template<> template<> struct OpInfo<IDN>::SiteOp<HAM>           { const static OpType Type = HAM; };
template<> template<> struct OpInfo<IDN>::SiteOp<CRE>           { const static OpType Type = CRE; };
template<> template<> struct OpInfo<IDN>::SiteOp<DES>           { const static OpType Type = DES; };
template<> template<> struct OpInfo<IDN>::SiteOp<CRE_COMP>      { const static OpType Type = CRE_COMP; };
template<> template<> struct OpInfo<IDN>::SiteOp<DES_COMP>      { const static OpType Type = DES_COMP; };

template<> template<> struct OpInfo<HAM>::SiteOp<HAM>           { const static OpType Type = IDN; };

template<> template<> struct OpInfo<CRE>::SiteOp<HAM>           { const static OpType Type = DES_COMP; };
template<> template<> struct OpInfo<CRE>::SiteOp<CRE>           { const static OpType Type = IDN; };

template<> template<> struct OpInfo<DES>::SiteOp<HAM>           { const static OpType Type = CRE_COMP; };
template<> template<> struct OpInfo<DES>::SiteOp<DES>           { const static OpType Type = IDN; };

template<> template<> struct OpInfo<CRE_COMP>::SiteOp<HAM>      { const static OpType Type = DES; };
template<> template<> struct OpInfo<CRE_COMP>::SiteOp<CRE_COMP> { const static OpType Type = IDN; };

template<> template<> struct OpInfo<DES_COMP>::SiteOp<HAM>      { const static OpType Type = CRE; };
template<> template<> struct OpInfo<DES_COMP>::SiteOp<DES_COMP> { const static OpType Type = IDN; };

std::string convert(OpType type) {
  switch(type & MASK) {
    case ZERO:
      return std::string("Null");
    case IDN:
      return std::string("Identity");
    case HAM:
      return std::string("Local Ham");
    case CRE:
      return std::string("Creation");
    case DES:
      return std::string("Destruction");
    case CRE_COMP:
      return std::string("Complement Creation");
    case DES_COMP:
      return std::string("Complement Destruction");
    default:
      return std::string("Unknown OpType");
  }
}

OpType get_site_optype(const OpType L, const OpType R) { return OpInfo<L>::SiteOp<R>::Type; }

int main() {
  OpType S_Op = get_site_optype(IDN, HAM); std::cout << convert(S_Op) << std::endl;
//OpType S_Op = OpInfo<IDN>::SiteOp<IDN>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<IDN>::SiteOp<HAM>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<IDN>::SiteOp<CRE>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<IDN>::SiteOp<DES>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<HAM>::SiteOp<IDN>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<HAM>::SiteOp<HAM>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<HAM>::SiteOp<CRE>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<HAM>::SiteOp<DES>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<CRE>::SiteOp<IDN>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<CRE>::SiteOp<HAM>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<CRE>::SiteOp<CRE>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<CRE>::SiteOp<DES>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<DES>::SiteOp<IDN>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<DES>::SiteOp<HAM>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<DES>::SiteOp<CRE>::Type; std::cout << convert(S_Op) << std::endl;
//       S_Op = OpInfo<DES>::SiteOp<DES>::Type; std::cout << convert(S_Op) << std::endl;
  return 0;
}
