#ifndef _MPSXX_CXX11_MP_STATES_H
#define _MPSXX_CXX11_MP_STATES_H 1

#include <vector>
#include <string>
#include <sstream>

#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx     {

enum MPS_TYPE { WAVEFUNCTION, LEFTCANONICAL, RIGHTCANONICAL };

template<class Q>
using MpStates = std::vector<btas::QSDArray<3, Q>>;

}; // namespace mpsxx

#endif // _MPSXX_CXX11_MP_STATES_H
