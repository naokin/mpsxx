
#include <sstream>

#include "bit_operator_type.h"

namespace mpsxx     {

namespace fermionic {

namespace bit_operator {
//                                                                                    + complementary
//                                                                                    |+ I or H
//                                                                                    ||   + 1st type      + 2nd type
//                                                                                    ||   |  + 1st index  |  + 2nd index
//                                                                                    ||   |  |            |  |
const BIT_OPERATOR_TYPE  ZERO        = 0x00000000; // 0000 00 000000000000 00 000000000000 -- zero component

// Masks

const BIT_OPERATOR_TYPE  IDEN        = 0x40000000; // 0100 00 000000000000 00 000000000000 -- identity operator
const BIT_OPERATOR_TYPE  HAM         = 0xc0000000; // 1100 00 000000000000 00 000000000000 -- hamiltonian

const BIT_OPERATOR_TYPE  TYPE        = 0xf0000000; // 1111 00 000000000000 00 000000000000 -- operator type mask
const BIT_OPERATOR_TYPE  NORMAL      = 0x7fffffff; // 0111 11 111111111111 11 111111111111 -- convert to normal operator
const BIT_OPERATOR_TYPE  COMP        = 0x80000000; // 1000 00 000000000000 00 000000000000 -- complementary operator mask

const BIT_OPERATOR_TYPE  MASK        = 0x7c003000; // 0111 11 000000000000 11 000000000000 -- operator mask
const BIT_OPERATOR_TYPE  FIRST       = 0x0fffc000; // 0000 11 111111111111 00 000000000000 -- first operator mask
const BIT_OPERATOR_TYPE  SECOND      = 0x00003fff; // 0000 00 000000000000 11 111111111111 -- second operator mask
const BIT_OPERATOR_TYPE  INDEX       = 0x03ffcfff; // 0000 00 111111111111 00 111111111111 -- index mask

const BIT_OPERATOR_TYPE  DOUBLE      = 0x30000000; // 0011 00 000000000000 00 000000000000 -- two index operator mask
const BIT_OPERATOR_TYPE  SINGLE      = 0x20000000; // 0010 00 000000000000 00 000000000000 -- one index operator mask 1
const BIT_OPERATOR_TYPE  SINGLE_2    = 0x10000000; // 0001 00 000000000000 00 000000000000 -- one index operator mask 2

const BIT_OPERATOR_TYPE  CONJ_D      = 0x08002000; // 0000 10 000000000000 10 000000000000 -- with XOR, conjugate two index operator
const BIT_OPERATOR_TYPE  CONJ_S_1    = 0x08000000; // 0000 10 000000000000 00 000000000000 -- with XOR, conjugate one index operator 1
const BIT_OPERATOR_TYPE  CONJ_S_2    = 0x00002000; // 0000 00 000000000000 10 000000000000 -- with XOR, conjugate one index operator 2

// Operator mask without spins

const BIT_OPERATOR_TYPE  CRE         = 0x28000000; // 0010 10 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  DES         = 0x20000000; // 0010 00 000000000000 00 000000000000

const BIT_OPERATOR_TYPE  CRE_CRE     = 0x38002000; // 0011 10 000000000000 10 000000000000
const BIT_OPERATOR_TYPE  CRE_DES     = 0x38000000; // 0011 10 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  DES_CRE     = 0x30002000; // 0011 00 000000000000 10 000000000000
const BIT_OPERATOR_TYPE  DES_DES     = 0x30000000; // 0011 00 000000000000 00 000000000000

// Operator mask with spins

const BIT_OPERATOR_TYPE  CRE_A_1     = 0x2c000000; // 0010 11 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  CRE_A_2     = 0x10003000; // 0001 00 000000000000 11 000000000000
const BIT_OPERATOR_TYPE  CRE_B_1     = 0x28000000; // 0010 10 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  CRE_B_2     = 0x10002000; // 0001 00 000000000000 10 000000000000

const BIT_OPERATOR_TYPE  DES_A_1     = 0x24000000; // 0010 01 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  DES_A_2     = 0x10001000; // 0001 00 000000000000 01 000000000000
const BIT_OPERATOR_TYPE  DES_B_1     = 0x20000000; // 0010 00 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  DES_B_2     = 0x10000000; // 0001 00 000000000000 00 000000000000

const BIT_OPERATOR_TYPE  CRE_A_CRE_A = 0x3c003000; // 0011 11 000000000000 11 000000000000
const BIT_OPERATOR_TYPE  CRE_A_CRE_B = 0x3c002000; // 0011 11 000000000000 10 000000000000
const BIT_OPERATOR_TYPE  CRE_B_CRE_A = 0x38003000; // 0011 10 000000000000 11 000000000000
const BIT_OPERATOR_TYPE  CRE_B_CRE_B = 0x38002000; // 0011 10 000000000000 10 000000000000

const BIT_OPERATOR_TYPE  CRE_A_DES_A = 0x3c001000; // 0011 11 000000000000 01 000000000000
const BIT_OPERATOR_TYPE  CRE_A_DES_B = 0x3c000000; // 0011 11 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  CRE_B_DES_A = 0x38001000; // 0011 10 000000000000 01 000000000000
const BIT_OPERATOR_TYPE  CRE_B_DES_B = 0x38000000; // 0011 10 000000000000 00 000000000000

const BIT_OPERATOR_TYPE  DES_A_CRE_A = 0x34003000; // 0011 01 000000000000 11 000000000000
const BIT_OPERATOR_TYPE  DES_A_CRE_B = 0x34002000; // 0011 01 000000000000 10 000000000000
const BIT_OPERATOR_TYPE  DES_B_CRE_A = 0x30003000; // 0011 00 000000000000 11 000000000000
const BIT_OPERATOR_TYPE  DES_B_CRE_B = 0x30002000; // 0011 00 000000000000 10 000000000000

const BIT_OPERATOR_TYPE  DES_A_DES_A = 0x34001000; // 0011 01 000000000000 01 000000000000
const BIT_OPERATOR_TYPE  DES_A_DES_B = 0x34000000; // 0011 01 000000000000 00 000000000000
const BIT_OPERATOR_TYPE  DES_B_DES_A = 0x30001000; // 0011 00 000000000000 01 000000000000
const BIT_OPERATOR_TYPE  DES_B_DES_B = 0x30000000; // 0011 00 000000000000 00 000000000000

const BIT_OPERATOR_TYPE  FIELD_SHIFT =   12;
const BIT_OPERATOR_TYPE  INDEX_SHIFT =   14;
const BIT_OPERATOR_TYPE  INDEX_LIMIT = 4096;

}; // namespace bit_operator
}; // namespace fermionic
}; // namespace mpsxx

mpsxx::fermionic::Quantum mpsxx::fermionic::get_quantum
(const mpsxx::fermionic::BIT_OPERATOR_TYPE& _op_type)
{
  using namespace bit_operator;
  int p = 0; // particle #
  int s = 0; // spin #
  if((_op_type & IDEN) == 0) { // except for I and H
    if(_op_type & SINGLE) {
      BIT_OPERATOR_TYPE op1 = (_op_type & MASK & FIRST)  >> INDEX_SHIFT >> FIELD_SHIFT;
      if(op1 & 2) { ++p; if(op1 & 1) ++s; else --s; }
      else        { --p; if(op1 & 1) --s; else ++s; }
    }
    if(_op_type & SINGLE_2) {
      BIT_OPERATOR_TYPE op2 = (_op_type & MASK & SECOND) >> FIELD_SHIFT;
      if(op2 & 2) { ++p; if(op2 & 1) ++s; else --s; }
      else        { --p; if(op2 & 1) --s; else ++s; }
    }
  }
  return std::move(Quantum(p, s));
}

std::string mpsxx::fermionic::translate
(const mpsxx::fermionic::BIT_OPERATOR_TYPE& _op_type)
{
  using namespace bit_operator;
  std::ostringstream opname;
  if(_op_type & IDEN) {
    if(_op_type & COMP) opname << "H";
    else             opname << "I";
  }
  else {
    if(_op_type & COMP) opname << "R:";
    if(_op_type & SINGLE) {
      BIT_OPERATOR_TYPE op1 = (_op_type & MASK & FIRST)  >> INDEX_SHIFT >> FIELD_SHIFT;
      if(op1 & 2) opname << "Cre";
      else        opname << "Des";
      if(op1 & 1) opname << "A[";
      else        opname << "B[";
      opname << ((_op_type & INDEX & FIRST) >> INDEX_SHIFT) << "]";
    }
    if(_op_type & SINGLE_2) {
      BIT_OPERATOR_TYPE op2 = (_op_type & MASK & SECOND)  >> FIELD_SHIFT;
      if(op2 & 2) opname << "_Cre";
      else        opname << "_Des";
      if(op2 & 1) opname << "A[";
      else        opname << "B[";
      opname << (_op_type & INDEX & SECOND) << "]";
    }
  }
  return opname.str();
}
