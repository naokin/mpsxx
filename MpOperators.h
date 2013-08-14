#ifndef _MPSXX_CXX11_MP_OPERATORS_H
#define _MPSXX_CXX11_MP_OPERATORS_H 1

#include <vector>
#include <string>
#include <sstream>

#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx     {

enum MPO_TYPE { HEISENBERG, FERMION_HUBBARD, MOLECULAR };

template<class Q>
using MpOperators = std::vector<btas::QSDArray<4, Q>>;

}; // namespace mpsxx

#endif // _MPSXX_CXX11_MP_OPERATORS_H
