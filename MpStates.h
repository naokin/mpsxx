#ifndef _MPSXX_CXX11_MP_STATES_H
#define _MPSXX_CXX11_MP_STATES_H 1

#include <vector>
#include <string>
#include <sstream>

#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx     {

template<class Q>
using MpStates = std::vector<btas::QSDArray<3, Q>>;

inline std::string get_mpsfile(const std::string& prefix, const int& index)
{
  std::stringstream filename;
  filename << prefix << "mps_site-" << index << /* mpigetrank() << */ ".tmp";
  return filename.str();
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_MP_STATES_H
