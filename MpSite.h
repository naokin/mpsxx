#ifndef __MPSXX_MPSITE_H
#define __MPSXX_MPSITE_H

#include <symmetry/fermion.h>

namespace mpsxx {

/// the struct MpSite determines the physical indices of every site
template<class Q>
struct MpSite {

  /// collection of physical indices
  enum PHYSICAL_INDEX { vacuum = 0 };

  /// \return collection of quantum numbers corresponding to the physical indices
  const static btas::Qshapes<Q> quanta() { return btas::Qshapes<Q>(1, Q::zero()); }

};

template<>
struct MpSite<fermion> {

  enum PHYSICAL_INDEX { vacuum = 0, alpha = 1, beta = 2, pair = 3 };

  /// \return the collection of physical quantumnumbers corresponding to the physical indices on a site.
  const static btas::Qshapes<fermion> quanta()
  {
    btas::Qshapes<fermion> q(4);
    q[vacuum] = fermion(0, 0);
    q[alpha ] = fermion(1,+1);
    q[beta  ] = fermion(1,-1);
    q[pair  ] = fermion(2, 0);
    return q;
  }

};

} // namespace mpsxx

#endif // __MPSXX_MPSITE_H
