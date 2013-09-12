#include <iostream>
#include <vector>
#include <cstring>

#include <btas/DENSE/DArray.h>

#include "generate_qc_operators.h"
#include "driver/parsing_integral.h"

int main(int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  using namespace mpsxx;

  bool enable_swap_sweep_dir = false;
  int iprint = 0;

  std::string f_dump = "FCIDUMP";
  std::string prefix = "./";
  std::string f_rord;

  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) f_dump = argv[++iarg];
    if(strcmp(argv[iarg],"-r") == 0) f_rord = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefix = argv[++iarg];
    if(strcmp(argv[iarg],"-x") == 0) enable_swap_sweep_dir = true;
    if(strcmp(argv[iarg],"-v") == 0) iprint = 1;
  }

  int Norbs;
  int Nelec;
  double Ecore;
  btas::DArray<2> oneint;
  btas::DArray<4> twoint;

  cout << "\t****************************************************************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::MPO GENERATOR FOR QC "                                                        << endl;
  cout << "\t****************************************************************************************************" << endl;
  cout << endl;
  cout << "\t====================================================================================================" << endl;
  cout << "\t\tLoading integral information from FCIDUMP "                                                         << endl;
  cout << "\t====================================================================================================" << endl;
  cout << endl;

  std::ifstream ist_dump(f_dump.c_str());
  if(!f_rord.empty()) {
    std::ifstream ist_rord(f_rord.c_str());
    std::vector<int> reorder;
    parsing_reorder(ist_rord, reorder);
    parsing_fcidump(ist_dump, Norbs, Nelec, Ecore, oneint, twoint, reorder);
  }
  else {
    parsing_fcidump(ist_dump, Norbs, Nelec, Ecore, oneint, twoint);
  }

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tGenerating QC MPOs from 1- and 2-particle integrals"                                                << endl;
  cout << "\t====================================================================================================" << endl;
  cout << endl;

  MpOperators<fermionic::Quantum> mpos(Norbs);
  fermionic::generate_qc_operators(mpos, oneint, twoint, enable_swap_sweep_dir, prefix);

  cout << endl;

  return 0;
}
