#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "make_sweep.h"
#include "make_random_mpss.h"
#include "gaugefix.hpp"

namespace mpsxx {

/// The main dmrg routine.
/// \param input DMRGInput object which contains the diferent parameters that define this dmrg run
std::vector<double> dmrg(const DMRGInput& input)
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  const size_t K = input.N_roots;
  const size_t MAX_ITER = input.N_max_sweep_iter;

  std::vector<double> esav(K,0.0);

  for(size_t iroot = 0; iroot < K; ++iroot) {
    if(!input.restart) make_random_mpss(input,iroot);
    bool conv = false;
    // Optimization with sweep algorithm
    for(size_t iter = 0; iter < MAX_ITER && !conv; ++iter) {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << endl;
      cout << "\t====================================================================================================" << endl;

      double eswp = make_sweep(esav,input,iroot);
      double edif = eswp-esav[iroot];

      cout << "\t====================================================================================================" << endl;
      cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << " FINISHED" << endl;
      cout.precision(16);
      cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
      cout.precision(2);
      cout << setw(8) << scientific << edif << " ) " << endl;
      cout << endl;

      esav[iroot] = eswp;
      if(iter > 0 && std::fabs(edif) < input.tolerance) conv = true;
    }
    // Stop by no convergence
    if(!conv) {
      cout << "\t+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-" << endl;
      cout << "\t\tNO CONVERGENCE MET FOR ROOT = " << setw(2) << iroot << endl;
      cout << "\t\tPROGRAM STOPPED..." << endl;
      break;
    }
    if(input.algorithm == ONESITE) continue;

    DMRGInput in2nd = input;
    in2nd.algorithm = ONESITE;
    in2nd.restart   = true;

    conv = false;
    // Optimization with sweep algorithm
    for(size_t iter = 0; iter < MAX_ITER && !conv; ++iter) {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << endl;
      cout << "\t====================================================================================================" << endl;

      double eswp = make_sweep(esav,in2nd,iroot);
      double edif = eswp-esav[iroot];

      cout << "\t====================================================================================================" << endl;
      cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << " FINISHED" << endl;
      cout.precision(16);
      cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
      cout.precision(2);
      cout << setw(8) << scientific << edif << " ) " << endl;
      cout << endl;

      esav[iroot] = eswp;
      if(iter > 0 && std::fabs(edif) < in2nd.tolerance) conv = true;
    }
    // Stop by no convergence
    if(!conv) {
      cout << "\t+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-" << endl;
      cout << "\t\tNO CONVERGENCE MET FOR ROOT = " << setw(2) << iroot << endl;
      cout << "\t\tPROGRAM STOPPED..." << endl;
      break;
    }
    gaugefix<fermion>(input,iroot);
  }
      cout << "\t====================================================================================================" << endl;
      cout.precision(16);
  for(size_t iroot = 0; iroot < K; ++iroot)
      cout << "\t\t\tSweep Energy = " << setw(24) << fixed << esav[iroot] << endl;

  return esav;
}

} // namespace mpsxx

int main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;

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
