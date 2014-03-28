#ifndef PROTOTYPE_MPSITE_H
#define PROTOTYPE_MPSITE_H

#include <vector>
#include <btas/DArray.h>
#include <btas/Dblas.h>

//
// workspaces of dmrg optimization
//
struct MpSite
{
  // hamiltonian
  btas::DArray<4>
    mpopr;
  // wavefunction and rotation matrices
  // used as bra states
  btas::DArray<3>
    wfncn;
  btas::DArray<3>
    lmpsn;
  btas::DArray<3>
    rmpsn;
  // renormalized operators
  btas::DArray<3>
    lstrn;
  btas::DArray<3>
    rstrn;

  // wavefunction and rotation matrices from other calculation
  // used as ket states
  std::vector< btas::DArray<3> >
    wfnc0;
  std::vector< btas::DArray<3> >
    lmps0;
  std::vector< btas::DArray<3> >
    rmps0;

  // renormalized overlaps and operators between current state and other states
  std::vector< btas::DArray<2> >
    lovl0n;
  std::vector< btas::DArray<2> >
    rovl0n;
  std::vector< btas::DArray<3> >
    lstr0n;
  std::vector< btas::DArray<3> >
    rstr0n;

  MpSite(int nroot) : wfnc0 (nroot, btas::DArray<3>()),
                      lmps0 (nroot, btas::DArray<3>()), rmps0 (nroot, btas::DArray<3>()),
                      lovl0n(nroot, btas::DArray<2>()), rovl0n(nroot, btas::DArray<2>()),
                      lstr0n(nroot, btas::DArray<3>()), rstr0n(nroot, btas::DArray<3>())
  {
  }

  void save(int iroot = 0)
  {
    btas::Dcopy(wfncn, wfnc0[iroot]);
    btas::Dcopy(lmpsn, lmps0[iroot]);
    btas::Dcopy(rmpsn, rmps0[iroot]);
  }
};

typedef std::vector<MpSite> MpStorages;

#endif // PROTOTYPE_MPSITE_H
