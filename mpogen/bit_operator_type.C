
#include <sstream>

#include "bit_operator_type.h"

fermion mpsxx::mpogen::get_quantum (const mpsxx::mpogen::BIT_OPERATOR_TYPE& _op_type)
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
  return std::move(fermion(p, s));
}

std::string mpsxx::mpogen::translate (const mpsxx::mpogen::BIT_OPERATOR_TYPE& _op_type)
{
  using namespace bit_operator;
  std::ostringstream opname;
  if(_op_type & IDEN) {
    if(_op_type & COMP) opname << "H";
    else             opname << "I";
  }
  else {
    if(_op_type & COMP) opname << "Comp::";
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
