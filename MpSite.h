#ifndef _MPSXX_CXX11_MPSITE_H
#define _MPSXX_CXX11_MPSITE_H 1

#include <symmetry/Fermion/Quantum.h>

namespace mpsxx {

   /**
    * the struct MpSite determines the physical indices of every site
    */
   template<class Q>
      struct MpSite {

         //!collection of physical indices
         enum PHYSICAL_INDEX { vacuum = 0 };

         //! @return collection of quantum numbers corresponding to the physical indices
         const static btas::Qshapes<Q> quanta() { return btas::Qshapes<Q>(1, Q::zero()); }

      };

   //! template Specializion for fermionic systems: 4 indices
   namespace fermionic { enum PHYSICAL_INDEX { vacuum = 0, alpha = 1, beta = 2, pair = 3 }; };

   template<>
      struct MpSite<fermionic::Quantum> {

         //!@return the collection of physical quantumnumbers corresponding to the physical indices on a site.
         const static btas::Qshapes<fermionic::Quantum> quanta() {
            btas::Qshapes<fermionic::Quantum> q(4);
            q[fermionic::vacuum] = fermionic::Quantum(0, 0);
            q[fermionic::alpha ] = fermionic::Quantum(1,+1);
            q[fermionic::beta  ] = fermionic::Quantum(1,-1);
            q[fermionic::pair  ] = fermionic::Quantum(2, 0);

            return q;

         }
      };

}; // namespace mpsxx

#endif // _MPSXX_CXX11_MPSITE_H
