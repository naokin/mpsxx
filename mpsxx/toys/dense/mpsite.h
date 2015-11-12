#ifndef _PROTOTYPE_MPSITE_H
#define _PROTOTYPE_MPSITE_H 1

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>
#include <legacy/DArray.h>

namespace prototype
{

struct MpSite
{
  //
  // storage label
  //
  std::string label;

  //
  // storage size
  //
  int nstore;

  //
  // matrix product operator (MPO)
  //
  btas::DArray<4> mpo;

  //
  // matrix product state (MPS) / left-canonical / right-canonical / wavefunction at this site
  //
  std::vector< btas::DArray<3> > lmps;
  std::vector< btas::DArray<3> > rmps;
  std::vector< btas::DArray<3> > wfnc;

  std::vector< btas::DArray<3> > ltnj;
  std::vector< btas::DArray<3> > rtnj;

  //
  // renormalized operator
  //
  std::vector< btas::DArray<3> > lopr;
  std::vector< btas::DArray<3> > ropr;

  //
  // file prefix
  //
  std::string fprefix;

  //
  // constructor
  //
  MpSite(const std::string& name = "mpsite", int nroots = 1, const std::string& dir = ".")
  : label(name),
    lmps(nroots, btas::DArray<3>()),
    rmps(nroots, btas::DArray<3>()),
    wfnc(nroots, btas::DArray<3>()),
    ltnj(nroots, btas::DArray<3>()),
    rtnj(nroots, btas::DArray<3>()),
    lopr(nroots, btas::DArray<3>()),
    ropr(nroots, btas::DArray<3>()),
    fprefix(dir)
  {
    nstore = nroots;
  }

  //
  // File-I/O
  //
  void load(int suffix)
  {
    std::ostringstream fname;
    fname << fprefix << "/" << label << "-" << suffix << ".tmp";
    std::ifstream fload(fname.str().c_str());
    boost::archive::binary_iarchive iarc(fload);

    iarc & nstore;
    iarc & mpo;

    iarc & lmps;
    iarc & rmps;
    iarc & wfnc;
    iarc & ltnj;
    iarc & rtnj;
    iarc & lopr;
    iarc & ropr;

    fload.close();
  }

  void save(int suffix)
  {
    std::ostringstream fname;
    fname << fprefix << "/" << label << "-" << suffix << ".tmp";
    std::ofstream fsave(fname.str().c_str());
    boost::archive::binary_oarchive oarc(fsave);

    oarc & nstore;
    oarc & mpo; mpo.free();

    oarc & lmps; lmps.clear();
    oarc & rmps; rmps.clear();
    oarc & wfnc; wfnc.clear();
    oarc & ltnj; ltnj.clear();
    oarc & rtnj; rtnj.clear();
    oarc & lopr; lopr.clear();
    oarc & ropr; ropr.clear();

    fsave.close();
  }
};

typedef std::vector<MpSite> MpStorages;

};

#endif // _PROTOTYPE_MPSITE_H
