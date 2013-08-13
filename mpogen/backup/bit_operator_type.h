#ifndef _MPSXX_CXX11_BIT_OPERATOR_TYPE_H
#define _MPSXX_CXX11_BIT_OPERATOR_TYPE_H 1

#include <cstring>

#include "Fermion/Quantum.h"

namespace mpsxx     {

namespace fermionic {

typedef unsigned int BIT_OPERATOR_TYPE;

namespace bit_operator {

const extern BIT_OPERATOR_TYPE  ZERO       ;

// Masks

const extern BIT_OPERATOR_TYPE  IDEN       ;
const extern BIT_OPERATOR_TYPE  HAM        ;

const extern BIT_OPERATOR_TYPE  TYPE       ;
const extern BIT_OPERATOR_TYPE  NORMAL     ;
const extern BIT_OPERATOR_TYPE  COMP       ;

const extern BIT_OPERATOR_TYPE  MASK       ;
const extern BIT_OPERATOR_TYPE  FIRST      ;
const extern BIT_OPERATOR_TYPE  SECOND     ;
const extern BIT_OPERATOR_TYPE  INDEX      ;

const extern BIT_OPERATOR_TYPE  DOUBLE     ;
const extern BIT_OPERATOR_TYPE  SINGLE     ;
const extern BIT_OPERATOR_TYPE  SINGLE_2   ;

const extern BIT_OPERATOR_TYPE  CONJ_D     ;
const extern BIT_OPERATOR_TYPE  CONJ_S_1   ;
const extern BIT_OPERATOR_TYPE  CONJ_S_2   ;

// Operator mask without spins

const extern BIT_OPERATOR_TYPE  CRE_1      ;
const extern BIT_OPERATOR_TYPE  CRE_2      ;
const extern BIT_OPERATOR_TYPE  DES_1      ;
const extern BIT_OPERATOR_TYPE  DES_2      ;

const extern BIT_OPERATOR_TYPE  CRE_CRE    ;
const extern BIT_OPERATOR_TYPE  CRE_DES    ;
const extern BIT_OPERATOR_TYPE  DES_CRE    ;
const extern BIT_OPERATOR_TYPE  DES_DES    ;

// Operator mask with spins

const extern BIT_OPERATOR_TYPE  CRE_A_1    ;
const extern BIT_OPERATOR_TYPE  CRE_A_2    ;
const extern BIT_OPERATOR_TYPE  CRE_B_1    ;
const extern BIT_OPERATOR_TYPE  CRE_B_2    ;

const extern BIT_OPERATOR_TYPE  DES_A_1    ;
const extern BIT_OPERATOR_TYPE  DES_A_2    ;
const extern BIT_OPERATOR_TYPE  DES_B_1    ;
const extern BIT_OPERATOR_TYPE  DES_B_2    ;

const extern BIT_OPERATOR_TYPE  CRE_A_CRE_A;
const extern BIT_OPERATOR_TYPE  CRE_A_CRE_B;
const extern BIT_OPERATOR_TYPE  CRE_B_CRE_A;
const extern BIT_OPERATOR_TYPE  CRE_B_CRE_B;

const extern BIT_OPERATOR_TYPE  CRE_A_DES_A;
const extern BIT_OPERATOR_TYPE  CRE_A_DES_B;
const extern BIT_OPERATOR_TYPE  CRE_B_DES_A;
const extern BIT_OPERATOR_TYPE  CRE_B_DES_B;

const extern BIT_OPERATOR_TYPE  DES_A_CRE_A;
const extern BIT_OPERATOR_TYPE  DES_A_CRE_B;
const extern BIT_OPERATOR_TYPE  DES_B_CRE_A;
const extern BIT_OPERATOR_TYPE  DES_B_CRE_B;

const extern BIT_OPERATOR_TYPE  DES_A_DES_A;
const extern BIT_OPERATOR_TYPE  DES_A_DES_B;
const extern BIT_OPERATOR_TYPE  DES_B_DES_A;
const extern BIT_OPERATOR_TYPE  DES_B_DES_B;

const extern BIT_OPERATOR_TYPE  FIELD_SHIFT;
const extern BIT_OPERATOR_TYPE  INDEX_SHIFT;
const extern BIT_OPERATOR_TYPE  INDEX_LIMIT;

}; // namespace bit_operator

Quantum get_quantum(const BIT_OPERATOR_TYPE& _type);
std::string translate(const BIT_OPERATOR_TYPE& _type);

}; // fermionic

}; // mpsxx

#endif // _MPSXX_CXX11_BIT_OPERATOR_TYPE_H
