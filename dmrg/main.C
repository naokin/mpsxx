#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
using namespace std;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; };

#include "input.h"
#include "mpsite.h"
#include "dmrg.h"
#include "lrt.h"
using namespace prototype;

int main(int argc, char* argv[])
{
  // by default
  DmrgInput dmrginput;

  int jobt = 1;

  string input;
  string output;
  string prefix;
  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) input  = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) output = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefix = argv[++iarg];
  }

  if(input.size() > 0) {
    ifstream fin(input.c_str());
    string entry;
    while(fin >> entry) {
      if(entry == "L")
        fin >> dmrginput.N_sites;
      if(entry == "nspin")
        fin >> dmrginput.N_spins;
      if(entry == "nelec")
        fin >> dmrginput.N_elecs;
      if(entry == "M")
        fin >> dmrginput.N_max_states;
      if(entry == "nroots")
        fin >> dmrginput.N_roots;
      if(entry == "tole")
        fin >> dmrginput.tolerance;
      if(entry == "Heisenberg" || entry == "heisenberg") {
        dmrginput.model = HEISENBERG;
        int    Nz = 1; // local Sz = Nz/2
        double J  = 1.0;
        double Jz = 1.0;
        double Hz = 0.0;
        while(fin >> entry) {
          if(entry == "Nz")
            fin >> Nz;
          else if(entry == "J" )
            fin >> J;
          else if(entry == "Jz")
            fin >> Jz;
          else if(entry == "Hz")
            fin >> Hz;
          else
            break;
        }
        dmrginput.heisenberg = HeisenbergModel(Nz, J, Jz, Hz);
      }
      if(entry == "Hubbard" || entry == "hubbard") {
        dmrginput.model = HUBBARD;
        double t = 1.0;
        double U = 1.0;
        while(fin >> entry) {
          if(entry == "t")
            fin >> t;
          else if(entry == "U")
            fin >> U;
          else
            break;
        }
        dmrginput.hubbard = HubbardModel(t, U);
      }
    }
  }
  if(prefix.size() > 0) {
    dmrginput.prefix = prefix;
  }

  //
  // assign cout as alias to fout
  //
  streambuf *backup;
  backup = cout.rdbuf();
  ofstream fout;
  if(output.size() > 0) {
    fout.open(output.c_str());
    cout.rdbuf(fout.rdbuf());
  }

  //
  // working spaces
  //

  MpStorages sites;

  //
  // dmrg optimization
  //

  cout << "\t****************************************************************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::OPTIMIZATION "                                                          << endl;
  cout << "\t****************************************************************************************************" << endl;

  dmrginput.energy = dmrg(dmrginput, sites);

//cout << "\t****************************************************************************************************" << endl;
//cout << "\t\t\t\tMPSXX::PROTOTYPE::TRANSFER OPERATOR ANALYSIS"                                                   << endl;
//cout << "\t****************************************************************************************************" << endl;

//analysis(dmrginput, sites);

  //
  // dmrg-lrt for excited states
  //

  if(dmrginput.N_roots > 1) {

    cout << "\t****************************************************************************************************" << endl;
    cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::LINEAR-RESPONSE "                                                       << endl;
    cout << "\t****************************************************************************************************" << endl;

    dmrg_lrt(dmrginput, sites);
  }

  cout << endl;

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
