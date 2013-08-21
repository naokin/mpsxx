
#include <time_stamp.h>

#include <symmetry/Fermion/Quantum.h>

#include "dmrg.h"
#include "initialize_mpstates.h"

#include "btas_template_specialize.h"

//
// random number generator
//
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX)-1.0; }

//
// DMRG main routine
//
double mpsxx::dmrg(const mpsxx::DmrgInput& input)
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::setprecision;
  using std::fixed;
  using std::scientific;

  const size_t N = input.N_sites;
  const int    M = input.N_max_states;

  MpOperators<fermionic::Quantum> mpos(N);
  MpStates   <fermionic::Quantum> mpss(N);

  if(!input.restart)
    initialize_mpstates(mpos, mpss, fermionic::Quantum(input.N_elecs, input.N_spins), rgen, input.prefix, input.N_max_states);

  double esav = 1.0e8;

  time_stamp ts;

  //
  // two-site optimization
  //
  for(size_t iter = 0; iter < 20; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(mpos, mpss, TWOSITE, input);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t----------------------------------------------------------------------------------------------------" << endl;
    cout << "\t\tTotal elapsed time: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. "          << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(std::fabs(edif) < input.tolerance) break;
  }

  //
  // one-site optimization
  //
  for(size_t iter = 0; iter < 20; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(mpos, mpss, ONESITE, input);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t----------------------------------------------------------------------------------------------------" << endl;
    cout << "\t\tTotal elapsed time: " << fixed << setprecision(2) << setw(8) << ts.elapsed() << " sec. "          << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(std::fabs(edif) < input.tolerance) break;
  }

  return esav;
}

