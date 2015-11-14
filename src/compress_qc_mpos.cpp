#include <iostream>
#include <iomanip>

#include "fileio.h"
#include "compress_qc_mpos.h"

size_t mpsxx::compress_mpos_2site (bool forward, btas::QSTArray<double,4,fermion>& x, btas::QSTArray<double,4,fermion>& y)
{
  btas::QSTArray<double,6,fermion> xy;
  btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,x,y,1.0,xy);

  btas::STArray<double,1> s;
  if(forward)
    btas::Gesvd<double,6,4,fermion,btas::LeftArrow> (xy,s,x,y,-6);
  else
    btas::Gesvd<double,6,4,fermion,btas::RightArrow>(xy,s,x,y,-6);

  size_t D = 0;
  std::cout << "block = ";
  for(auto it = s.begin(); it != s.end(); ++it) {
    size_t Di = it->second->size();
    std::cout << Di << ":";
    D += Di;
  }
  std::cout << "NET(" << D << ")" << std::endl;

  if(forward)
    btas::Dimm(s,y);
  else
    btas::Dimm(x,s);

  return D;
}

std::vector<size_t> mpsxx::compress_mpos_cycle (std::vector<btas::QSTArray<double,4,fermion>>& mpos)
{
  std::vector<size_t> Dsav(mpos.size()-1,0);

  size_t iter = 0;
  bool conv = false;
  while(!conv) {
    std::cout << "\t====================================================================================================" << std::endl;
    std::cout << "\t\tSWEEP :: " << std::setw(3) << iter << std::endl;
    std::cout << "\t----------------------------------------------------------------------------------------------------" << std::endl;
    // forward sweep
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "\t\t\tFORWARD SWEEP" << std::endl;
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::vector<size_t> Dfwd(mpos.size()-1,0);
    for(size_t i = 0; i < mpos.size()-1; ++i) {
      std::cout << "\t\t\tSITE [ " << std::setw(3) << i << " ] :: ";
      Dfwd[i] = compress_mpos_2site(1,mpos[i],mpos[i+1]);
    }
    // backward sweep
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "\t\t\tBACKWARD SWEEP" << std::endl;
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::vector<size_t> Dbwd(mpos.size()-1,0);
    for(size_t i = mpos.size()-1; i > 0; --i) {
      std::cout << "\t\t\tSITE [ " << std::setw(3) << i << " ] :: ";
      Dbwd[i-1] = compress_mpos_2site(0,mpos[i-1],mpos[i]);
    }
    conv  = std::equal(Dfwd.begin(),Dfwd.end(),Dsav.begin());
    conv &= std::equal(Dbwd.begin(),Dbwd.end(),Dsav.begin());
    ++iter;
    Dsav = Dbwd;
  }
  return Dsav;
}

std::vector<int> mpsxx::compress_qc_mpos (
  const std::vector<int>& groups,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& comp)
{
  const size_t MAXBUF = 200;

  size_t ig;
  std::set<int> ng(groups.begin(),groups.end());

  comp.resize(ng.size());

  std::vector<std::vector<size_t>> Ds(ng.size());

  std::cout << "\t\tThere are " << ng.size() << " groups to be separatory compressed." << std::endl;

  ig = 0;
  for(auto it = ng.begin(); it != ng.end(); ++it, ++ig) {
    int label = *it;
    std::cout << "\t\tCompressing operator group " << label << std::endl;

    size_t t = 0;
    for(; t < groups.size() && groups[t] != label; ++t);

    size_t ngen = 1;
    size_t nbuf = 1;
    comp[ig].swap(mpos[t]); ++t;

    size_t norb = comp[ig].size();

    for(; t < groups.size(); ++t) {
      if(groups[t] == label) {

        ++ngen;

        {
          btas::QSTArray<double,4,fermion> temp(comp[ig][0]);
          comp[ig][0].clear();
          btas::IVector<3> traceIdx = {0,1,2};
          btas::QSTdsum(temp,mpos[t][0],traceIdx,comp[ig][0]);
        }

        for(size_t s = 1; s < norb-1; ++s) {
          btas::QSTArray<double,4,fermion> temp(comp[ig][s]);
          comp[ig][s].clear();
          btas::IVector<2> traceIdx = {1,2};
          btas::QSTdsum(temp,mpos[t][s],traceIdx,comp[ig][s]);
        }

        {
          btas::QSTArray<double,4,fermion> temp(comp[ig][norb-1]);
          comp[ig][norb-1].clear();
          btas::IVector<3> traceIdx = {1,2,3};
          btas::QSTdsum(temp,mpos[t][norb-1],traceIdx,comp[ig][norb-1]);
        }
        // deallocation
        std::vector<btas::QSTArray<double,4,fermion>>().swap(mpos[t]);

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
//  if(nbuf > 1)
    Ds[ig] = compress_mpos_cycle(comp[ig]);
    std::cout << "\t\t\t";
    for(size_t s = 0; s < Ds[ig].size(); ++s) std::cout << Ds[ig][s] << ":";
    std::cout << ":compressed by " << ngen << " terms." << std::endl;
  }

  std::cout << "\t\tCheck boundary size..." << std::endl;
  ig = 0;
  for(auto it = ng.begin(); it != ng.end(); ++it, ++ig) {
    std::cout << "\t\t\t";
    for(size_t s = 0; s < Ds[ig].size(); ++s) std::cout << Ds[ig][s] << ":";
    std::cout << ":group = " << *it << std::endl;
  }

  return std::vector<int>(ng.begin(),ng.end());
}

void mpsxx::compress_qc_mpos (const size_t& N, const std::string& opname, const std::string& prefix)
{
  std::vector<size_t> Msav(N-1,0);

  size_t iter = 0;
  bool conv = false;
  while(!conv) {
    std::cout << "\t====================================================================================================" << std::endl;
    std::cout << "\t\tMPO compression sweep :: " << iter << std::endl;
    btas::QSTArray<double,4,fermion> lmpo;
    btas::QSTArray<double,4,fermion> rmpo;

    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "\t\t\tFORWARD SWEEP" << std::endl;
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    // forward sweep
    std::vector<size_t> Mfwd(N-1,0);
    load(lmpo,getfile(opname,prefix,0));
    for(size_t i = 0; i < N-1; ++i) {
      std::cout << "\t====================================================================================================" << std::endl;
      std::cout << "\t\tSITE [ " << std::setw(3) << i << " ] " << std::endl;
      std::cout << "\t----------------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "\t\t\t";

      load(rmpo,getfile(opname,prefix,i+1));
      Mfwd[i] = compress_mpos_2site(1,lmpo,rmpo);
      save(lmpo,getfile(opname,prefix,i));
      lmpo = rmpo;
    }
    save(lmpo,getfile(opname,prefix,N-1));

    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "\t\t\tBACKWARD SWEEP" << std::endl;
    std::cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    // backward sweep
    std::vector<size_t> Mbwd(N-1,0);
    for(size_t i = N-1; i > 0; --i) {
      std::cout << "\t====================================================================================================" << std::endl;
      std::cout << "\t\tSITE [ " << std::setw(3) << i << " ] " << std::endl;
      std::cout << "\t----------------------------------------------------------------------------------------------------" << std::endl;
      std::cout << "\t\t\t";

      load(lmpo,getfile(opname,prefix,i-1));
      Mbwd[i-1] = compress_mpos_2site(0,lmpo,rmpo);
      save(rmpo,getfile(opname,prefix,i));
      rmpo = lmpo;
    }
    save(rmpo,getfile(opname,prefix,0));

    conv  = std::equal(Mfwd.begin(),Mfwd.end(),Msav.begin());
    conv &= std::equal(Mbwd.begin(),Mbwd.end(),Msav.begin());

    ++iter;

    Msav = Mbwd;
  }

  size_t Msum = 0;
  std::cout << "\t\tCheck boundary size..." << std::endl;
  std::cout << "\t\t\tM = 1:";
  for(size_t s = 0; s < Msav.size(); ++s) {
    Msum += Msav[s];
    std::cout << Msav[s] << ":";
  }
  std::cout << "1" << std::endl;
}
