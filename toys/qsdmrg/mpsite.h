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

#include <vector>
#include <btas/QSDArray.h>

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
  btas::QSDArray<4> mpo;

  //
  // matrix product state (MPS) / left-canonical / right-canonical / wavefunction at this site
  //
  std::vector< btas::QSDArray<3> > lmps;
  std::vector< btas::QSDArray<3> > rmps;
  std::vector< btas::QSDArray<3> > wfnc;

  //
  // renormalized operator
  //
  std::vector< btas::QSDArray<3> > lopr;
  std::vector< btas::QSDArray<3> > ropr;

  //
  // file prefix
  //
  std::string fprefix;

  //
  // constructor
  //
  MpSite(const std::string& name = "mpsite", int nroots = 1, const std::string& dir = ".")
  : label(name),
    lmps(nroots),
    rmps(nroots),
    wfnc(nroots),
    lopr(nroots),
    ropr(nroots),
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
//  cout << "\t\t\tloading site from file, " << fname.str() << endl;
    std::ifstream fload(fname.str().c_str());
    boost::archive::binary_iarchive iarc(fload);

    int nroots; iarc >> nroots;

    if(nroots > nstore) {
      nstore = nroots;
      lmps.resize(nstore);
      rmps.resize(nstore);
      wfnc.resize(nstore);
      lopr.resize(nstore);
      ropr.resize(nstore);
    } // otherwise retain nstore

    iarc >> mpo;

    for(int i = 0; i < nroots; ++i) {
      iarc >> lmps[i];
      iarc >> rmps[i];
      iarc >> wfnc[i];
      iarc >> lopr[i];
      iarc >> ropr[i];
    }
  }

  void save(int suffix)
  {
    std::ostringstream fname;
    fname << fprefix << "/" << label << "-" << suffix << ".tmp";
//  cout << "\t\t\tsaving site into file, " << fname.str() << endl;
    std::ofstream fsave(fname.str().c_str());
    boost::archive::binary_oarchive oarc(fsave);

    oarc << nstore;

    oarc << mpo; mpo.clear();

    for(int i = 0; i < nstore; ++i) {
      oarc << lmps[i]; lmps[i].clear();
      oarc << rmps[i]; rmps[i].clear();
      oarc << wfnc[i]; wfnc[i].clear();
      oarc << lopr[i]; lopr[i].clear();
      oarc << ropr[i]; ropr[i].clear();
    }
  }

  void free()
  {
    mpo.clear();
    for(int i = 0; i < nstore; ++i) {
      lmps[i].clear();
      rmps[i].clear();
      wfnc[i].clear();
      lopr[i].clear();
      ropr[i].clear();
    }
  }
};

typedef std::vector<MpSite> MpStorages;

};

#endif // _PROTOTYPE_MPSITE_H
