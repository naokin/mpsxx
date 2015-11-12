#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include <legacy/DENSE/TArray.h>

#include "fermion.h"
#include "expectation.hpp"

int main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;

  mpsxx::DMRGInput input;

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
      if(entry == "restart")
        input.restart = true;
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
      if(entry == "onesite" || entry == "onedot")
        input.algorithm = mpsxx::ONESITE;
      if(entry == "twosite" || entry == "twodot")
        input.algorithm = mpsxx::TWOSITE;
      if(entry == "maxiter")
        fin >> input.N_max_sweep_iter;
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
  cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::EXPECTATION "                                                           << endl;
  cout << "\t****************************************************************************************************" << endl;

  size_t K = input.N_roots;
  btas::TArray<double,2> H(K,K);
  btas::TArray<double,2> S(K,K);

  for(size_t i = 0; i < K; ++i) {
    H(i,i) = mpsxx::expectation<fermion>(input,i,i);
    S(i,i) = mpsxx::overlap<fermion>(input,i,i);
    for(size_t j = 0; j < i; ++j) {
      double Hij = mpsxx::expectation<fermion>(input,i,j);
      H(i,j) = Hij;
      H(j,i) = Hij;
      double Sij = mpsxx::overlap<fermion>(input,i,j);
      S(i,j) = Sij;
      S(j,i) = Sij;
    }
  }

  cout.precision(8);
  cout << "Effective Hamiltonian :: " << endl;
  for(size_t i = 0; i < K; ++i) {
    for(size_t j = 0; j < K; ++j) {
      cout << fixed << setw(12) << H(i,j);
    }
    cout << endl;
  }
  cout << "Overlap Matrix :: " << endl;
  for(size_t i = 0; i < K; ++i) {
    for(size_t j = 0; j < K; ++j) {
      cout << fixed << setw(12) << S(i,j);
    }
    cout << endl;
  }

  btas::TArray<double,1> E;
  btas::TArray<double,2> C;
  btas::Sygv(1,'V','U',H,S,E,C);
  cout << "Eigenvalues :: " << endl;
  for(size_t i = 0; i < K; ++i) {
    cout << fixed << setw(12) << E(i) << endl;
  }

  cout << endl;

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
