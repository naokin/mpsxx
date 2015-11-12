#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include <time_stamp.h>

#include "mpidefs.h" // Communicator

#include "make_sweep.h"
#include "make_random_mpss.h"
#include "gaugefix.hpp"

namespace mpsxx {

/// The main dmrg routine.
/// \param input DMRGInput object which contains the diferent parameters that define this dmrg run
std::vector<double> dmrg (const DMRGInput& input)
{
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  const size_t K_roots = input.N_roots;
  const size_t MAX_ITER = input.N_max_sweep_iter;

  std::vector<double> esav(K_roots,0.0);

  Communicator world;

  for(size_t iroot = 0; iroot < K_roots; ++iroot) {

    // Build initial MPSs
    if(!input.restart) make_random_mpss(input,iroot);

    bool conv = false;
    // Optimization with sweep algorithm
    for(size_t iter = 0; iter < MAX_ITER && !conv; ++iter) {

      pout << "\t====================================================================================================" << endl;
      pout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << endl;
      pout << "\t====================================================================================================" << endl;

      double eswp = make_sweep(esav,input,iroot);

      double edif = eswp-esav[iroot];
      pout << "\t====================================================================================================" << endl;
      pout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << " FINISHED" << endl;
      pout.precision(16);
      pout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
      pout.precision(2);
      pout << setw(8) << scientific << edif << " ) " << endl;
      pout << endl;
      if(world.rank() == 0) {
        esav[iroot] = eswp;
        if(iter > 0 && std::fabs(edif) < input.tolerance) conv = true;
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,conv,0);
      boost::mpi::broadcast(world,esav,0);
#endif
    }
    // Stop by no convergence
    if(!conv) {
      pout << "\t+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-" << endl;
      pout << "\t\tNO CONVERGENCE MET FOR ROOT = " << setw(2) << iroot << endl;
      pout << "\t\tPROGRAM STOPPED..." << endl;
      break;
    }

    if(input.algorithm == ONESITE) continue;

    DMRGInput in2nd = input;
    in2nd.algorithm = ONESITE;
    in2nd.restart   = true;

    conv = false;
    // Optimization with sweep algorithm
    for(size_t iter = 0; iter < MAX_ITER && !conv; ++iter) {

      pout << "\t====================================================================================================" << endl;
      pout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << endl;
      pout << "\t====================================================================================================" << endl;

      double eswp = make_sweep(esav,in2nd,iroot);

      double edif = eswp-esav[iroot];
      pout << "\t====================================================================================================" << endl;
      pout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] :: ROOT = " << setw(2) << iroot << " FINISHED" << endl;
      pout.precision(16);
      pout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
      pout.precision(2);
      pout << setw(8) << scientific << edif << " ) " << endl;
      pout << endl;
      if(world.rank() == 0) {
        esav[iroot] = eswp;
        if(iter > 0 && std::fabs(edif) < in2nd.tolerance) conv = true;
      }
#ifndef _SERIAL
      boost::mpi::broadcast(world,conv,0);
      boost::mpi::broadcast(world,esav,0);
#endif
    }
    // Stop by no convergence
    if(!conv) {
      pout << "\t+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-" << endl;
      pout << "\t\tNO CONVERGENCE MET FOR ROOT = " << setw(2) << iroot << endl;
      pout << "\t\tPROGRAM STOPPED..." << endl;
      break;
    }
//  gaugefix<fermion>(input,iroot);
  }
  pout << "\t====================================================================================================" << endl;
  pout.precision(16);
  for(size_t iroot = 0; iroot < K_roots; ++iroot) {
    pout << "\t\t\tSweep Energy = " << setw(24) << fixed << esav[iroot] << endl;
  }

  return esav;
}

} // namespace mpsxx

int main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;

#ifndef _SERIAL
  boost::mpi::environment env(argc,argv);
#endif
  Communicator world;

  std::string f_inp = "dmrg.conf";
  std::string f_out;
  std::string prefx = ".";

  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-i") == 0) f_inp = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) f_out = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefx = argv[++iarg];
  }

  mpsxx::DMRGInput input(f_inp); input.prefix = prefx;

  //
  // assign cout as alias to fout
  //
  std::streambuf *backup;
  backup = cout.rdbuf();
  std::ofstream fout;
  if(f_out.size() > 0) {
    std::ostringstream oss;
    oss << f_out << "." << world.rank();
    fout.open(oss.str().c_str());
    cout.rdbuf(fout.rdbuf());
  }

  time_stamp ts;

  //
  // dmrg optimization
  //
  input.energy = mpsxx::dmrg(input);

  pout << endl;
  pout.precision(2);
  pout << "\t\t\tTotal elapsed time: " << setw(8) << fixed << ts.elapsed() << endl;

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
