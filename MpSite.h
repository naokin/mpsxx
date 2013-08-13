#ifndef _MPSXX_CXX11_MPSITE_H
#define _MPSXX_CXX11_MPSITE_H 1

#include <symmetry/Fermion/Quantum.h>

namespace mpsxx {

template<class Q>
struct MpSite {
  enum PHYSICAL_INDEX { vacuum = 0 };
  const static btas::Qshapes<Q> quanta() { return btas::Qshapes<Q>(1, Q::zero()); }
};

template<>
struct MpSite<fermionic::Quantum> {
  enum PHYSICAL_INDEX { vacuum = 0, alpha = 1, beta = 2, pair = 3 };
  const static btas::Qshapes<fermionic::Quantum> quanta() {
    btas::Qshapes<fermionic::Quantum> q(4);
    q[vacuum] = fermionic::Quantum(0, 0);
    q[alpha ] = fermionic::Quantum(1,+1);
    q[beta  ] = fermionic::Quantum(1,-1);
    q[pair  ] = fermionic::Quantum(2, 0);
    return q;
  }
};

}; // namespace mpsxx

#endif // _MPSXX_CXX11_MPSITE_H
