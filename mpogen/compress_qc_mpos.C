#include <iostream>
#include <iomanip>

#include "compress_qc_mpos.h"

size_t mpsxx::compress_mpos_2site (bool forward, MPO<double,fermion>& x, MPO<double,fermion>& y)
{
  btas::QSTArray<double,6,fermion> xy;
  btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,x,y,1.0,xy);

  btas::STArray<double,1> s;
  btas::Gesvd(xy,s,x,y,-8);

  size_t D = 0;
  for(auto it = s.begin(); it != s.end(); ++it) D += it->second->size();

  if(forward)
    btas::Dimm(s,y);
  else
    btas::Dimm(x,s);

  return D;
}

std::vector<size_t> mpsxx::compress_mpos_cycle (std::vector<mpsxx::MPO<double,fermion>>& mpos)
{
  std::vector<size_t> Dfwd(mpos.size()-1,0);
  std::vector<size_t> Dbwd(mpos.size()-1,0);

  size_t iter = 0;
  bool conv = false;
  while(!conv) {
    std::cout << "\t\t\t\tsweep :: " << iter << std::endl;
    // forward sweep
    for(size_t i = 0; i < mpos.size()-1; ++i) {
      Dfwd[i] = compress_mpos_2site(1,mpos[i],mpos[i+1]);
    }
    // backward sweep
    for(size_t i = mpos.size()-1; i > 0; --i) {
      Dbwd[i-1] = compress_mpos_2site(0,mpos[i-1],mpos[i]);
    }
    conv = std::equal(Dfwd.begin(),Dfwd.end(),Dbwd.begin());
    ++iter;
  }
  return Dfwd;
}

void mpsxx::compress_qc_mpos (
  const std::vector<int>& groups,
  const std::vector<std::vector<mpsxx::MPO<double,fermion>>>& mpos,
        std::vector<std::vector<mpsxx::MPO<double,fermion>>>& comp)
{
  const size_t MAXBUF = 200;

  size_t ig;
  std::set<int> ng(groups.begin(),groups.end());

  comp.resize(ng.size());

  std::vector<std::vector<size_t>> Ds(ng.size());

  ig = 0;
  for(auto it = ng.begin(); it != ng.end(); ++it, ++ig) {
    int label = *it;
    std::cout << "\t\tCompressing operator group " << label << std::endl;

    size_t t = 0;
    for(; t < groups.size() && groups[t] != label; ++t);

    size_t nbuf = 1;
    comp[ig] = mpos[t]; ++t;

    size_t norb = comp[ig].size();

    for(; t < groups.size(); ++t) {
      if(groups[t] == label) {

        {
          MPO<double,fermion> temp(comp[ig][0]);
          comp[ig][0].clear();
          btas::IVector<3> traceIdx = {0,1,2};
          btas::QSTdsum(temp,mpos[t][0],traceIdx,comp[ig][0]);
        }

        for(size_t s = 1; s < norb-1; ++s) {
          MPO<double,fermion> temp(comp[ig][s]);
          comp[ig][s].clear();
          btas::IVector<2> traceIdx = {1,2};
          btas::QSTdsum(temp,mpos[t][s],traceIdx,comp[ig][s]);
        }

        {
          MPO<double,fermion> temp(comp[ig][norb-1]);
          comp[ig][norb-1].clear();
          btas::IVector<3> traceIdx = {1,2,3};
          btas::QSTdsum(temp,mpos[t][norb-1],traceIdx,comp[ig][norb-1]);
        }

        if(++nbuf >= MAXBUF) {
          std::cout << "\t\t\tdo compress mpos up to " << std::setw(6) << t << std::endl;
          Ds[ig] = compress_mpos_cycle(comp[ig]);
          std::cout << "\t\t\t";
          for(size_t s = 0; s < Ds[ig].size(); ++s) std::cout << Ds[ig][s] << ":";
          std::cout << std::endl;
          nbuf = 1;
        }
      }
    }
    if(nbuf > 1)
      Ds[ig] = compress_mpos_cycle(comp[ig]);
  }

  std::cout << "\t\tCheck boundary size..." << std::endl;
  ig = 0;
  for(auto it = ng.begin(); it != ng.end(); ++it, ++ig) {
    std::cout << "\t\t\t";
    for(size_t s = 0; s < Ds[ig].size(); ++s) std::cout << Ds[ig][s] << ":";
    std::cout << ":group = " << *it << std::endl;
  }
}
