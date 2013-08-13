#ifndef _MPSXX_CXX11_MP_OPERATORS_H
#define _MPSXX_CXX11_MP_OPERATORS_H 1

#include <vector>
#include <string>
#include <sstream>

#include <btas/QSPARSE/QSDArray.h>

namespace mpsxx     {

template<class Q>
using MpOperators = std::vector<btas::QSDArray<4, Q>>;

inline std::string get_mpofile(const std::string& prefix, const int& index)
{
  std::stringstream filename;
  filename << prefix << "mpo_site-" << index << /* mpigetrank() << */ ".tmp";
  return filename.str();
}

}; // namespace mpsxx

#endif // _MPSXX_CXX11_MP_OPERATORS_H
