#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include "dmrg.h"
using namespace std;
using namespace btas;

// random # generator
double random_gen()
{
  return (2.0*rand()-1.0)/RAND_MAX;
}

int main(int argc, char* argv[])
{
  // by default
  DmrgInput dmrginput;

  string input;
  string output;
  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) input   = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) output  = argv[++iarg];
  }
  double Sz_max = 0.5;
  if(input.size() > 0) {
    ifstream fin(input.c_str());
    string entry;
    while(fin >> entry) {
      if(entry == "L" ) fin >> dmrginput.N_sites;
      if(entry == "Sz") fin >> Sz_max;
      if(entry == "M" ) fin >> dmrginput.N_max_states;
      if(entry == "J" ) fin >> dmrginput.J;
      if(entry == "Jz") fin >> dmrginput.Jz;
      if(entry == "hz") fin >> dmrginput.hz;
    }
  }
  dmrginput.N_phys_index = static_cast<int>(Sz_max * 2 + 1.5);

  streambuf *backup;
  backup = cout.rdbuf();
  ofstream fout;
  if(output.size() > 0) {
    fout.open(output.c_str());
    cout.rdbuf(fout.rdbuf());
  }

  cout << dmrginput << endl;

  //
  // working spaces
  //

  int nroot = 3;
  int L = dmrginput.N_sites;
  MpStorages sites(L, MpSite(nroot));

  //
  // dmrg optimization
  //

  vector<double> energy;
  energy.reserve(nroot);
  for(int iroot = 0; iroot < nroot; ++iroot) {
    cout << "\t**************************************************"
         <<   "**************************************************" << endl;
    cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::OPTIMIZATION "
         <<         "ROOT [ " << setw(2) << iroot << " ]  "        << endl;
    cout << "\t**************************************************"
         <<   "**************************************************" << endl;

    cout << "\t\t\t\t\tROOT [ " << setw(2) << iroot << " ] "       << endl;
    cout << "\t**************************************************"
         <<   "**************************************************" << endl;

    energy.push_back(dmrg(cout, dmrginput, sites, energy));

    for(int i = 0; i < L; ++i) sites[i].save(iroot);
  }

  cout << endl;

  cout.precision(4);
  for(int i = 0; i < sites.size()-1; ++i) {
    for(int iroot = 0; iroot < nroot; ++iroot) {
      for(int jroot = iroot+1; jroot < nroot; ++jroot) {
        DArray<2> scr;
        Dgemm(Trans, NoTrans, 1.0, sites[i].lmps0[iroot], sites[i].lmps0[jroot], 1.0, scr);
        cout << "sites[" << setw(2) << i << "]:" << setw(2) << iroot << "-" << setw(2) << jroot << endl;
        cout << scr << endl;
      }
    }
  }

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
