#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>

#include <legacy/DENSE/TArray.h>
#include <legacy/QSPARSE/QSTArray.h>

#include "parsing_integral.h"
#include "fileio.h"

#include "gen_qc_naive_mpos.h"
#include "compress_qc_mpos.h"

int main(int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  using namespace mpsxx;

  std::string fcidmp = "FCIDUMP";
  std::string reordr;
  std::string prefix = ".";

  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) fcidmp = argv[++iarg];
    if(strcmp(argv[iarg],"-r") == 0) reordr = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefix = argv[++iarg];
  }

  int Norbs;
  int Nelec;
  double Ecore;
  btas::TArray<double,2> oneint;
  btas::TArray<double,4> twoint;

  cout << "\t****************************************************************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::MPO GENERATOR FOR QC "                                                        << endl;
  cout << "\t****************************************************************************************************" << endl;
  cout << endl;
  cout << "\t====================================================================================================" << endl;
  cout << "\t\tLoading integral information from FCIDUMP "                                                         << endl;
  cout << "\t====================================================================================================" << endl;
  cout << endl;

  std::ifstream ifdump(fcidmp.c_str());
  if(!reordr.empty()) {
    std::ifstream ireord(reordr.c_str());
    std::vector<int> reorder;
    parsing_reorder(ireord, reorder);
    parsing_fcidump(ifdump, Norbs, Nelec, Ecore, oneint, twoint, reorder);
  }
  else {
    parsing_fcidump(ifdump, Norbs, Nelec, Ecore, oneint, twoint);
  }

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tGenerating QC MPOs from 1- and 2-particle integrals"                                                << endl;
  cout << "\t====================================================================================================" << endl;
  cout << endl;

  std::vector<std::vector<btas::QSTArray<double,4,fermion>>> mpos;
  std::vector<int> groups = gen_qc_naive_mpos(Norbs, oneint, twoint, mpos);

  std::cout << "\t\t" << std::setw(6) << mpos.size() << " operators have generated." << std::endl;
  std::fill(groups.begin(),groups.end(),1);

  std::vector<std::vector<btas::QSTArray<double,4,fermion>>> comp;
  compress_qc_mpos(groups,mpos,comp);

  for(size_t i = 0; i < comp[0].size(); ++i) save(comp[0][i], getfile("mpo",prefix,i));

  return 0;
}
