#include <iostream>
#include <vector>
#include <cstring>

#include <btas/DENSE/DArray.h>

#include "generate_qc_operators.h"
#include "parsing_fcidump.h"

int main(int argc, char* argv[])
{
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

  MpOperators<fermionic::Quantum> mpos(Norbs);
  fermionic::generate_qc_operators(mpos, oneint, twoint, enable_swap_sweep_dir, prefix);

  return 0;
}
