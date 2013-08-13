#ifndef _MPSXX_MP_OPERATORS_H
#define _MPSXX_MP_OPERATORS_H 1

#include <vector>

#include <btas/QSPARSE/QSDArray.h>

#include "Fermionic/Quantum.h"

namespace mpsxx {

typedef std::vector<btas::QSDArray<4, fermionic::Quantum>> MpOperators;

};

#endif // _MPSXX_MP_OPERATOR_H
