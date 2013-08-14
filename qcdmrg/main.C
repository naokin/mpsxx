#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "dmrg.h"

int main(int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  mpsxx::DmrgInput input;

  std::string f_inp = "dmrg.conf";
  std::string f_out;

  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-i") == 0) f_inp = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) f_out = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) input.prefix = argv[++iarg];
  }

  if(f_inp.size() > 0) {
    std::ifstream fin(f_inp.c_str());
    std::string entry;
    while(fin >> entry) {
      if(entry == "N")
        fin >> input.N_sites;
      if(entry == "spin")
        fin >> input.N_spins;
      if(entry == "elec")
        fin >> input.N_elecs;
      if(entry == "M" || entry == "max_states")
        fin >> input.N_max_states;
      if(entry == "nroots")
        fin >> input.N_roots;
      if(entry == "tole" || entry == "tolerance")
        fin >> input.tolerance;
    }
  }

  //
  // assign cout as alias to fout
  //
  std::streambuf *backup;
  backup = cout.rdbuf();
  std::ofstream fout;
  if(f_out.size() > 0) {
    fout.open(f_out.c_str());
    cout.rdbuf(fout.rdbuf());
  }

  cout << "\t****************************************************************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::OPTIMIZATION "                                                          << endl;
  cout << "\t****************************************************************************************************" << endl;

  //
  // dmrg optimization
  //
  input.energy = mpsxx::dmrg(input);

  cout << endl;

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
