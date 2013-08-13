#ifndef _MPSXX_MP_OPERATORS_H
#define _MPSXX_MP_OPERATORS_H 1

#include <vector>

#include <btas/QSPARSE/QSDArray.h>

#include "Fermion/Quantum.h"

namespace mpsxx     {

namespace fermionic {

typedef std::vector<btas::QSDArray<4, Quantum>> MpOperators;

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_MP_OPERATOR_H
