#include <symmetry/Fermion/Quantum.h>

#include "dmrg.h"
#include "initialize_mpstates.h"

#include "btas_template_specialize.h"

/**
 * simple random number generator
*/
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX)-1.0; }

/**
 * The main dmrg routine.
 * @param input DmrgInput object which contains the diferent parameters that define this dmrg run
 */
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

  MPO<fermionic::Quantum> mpos(N);
  MPS<fermionic::Quantum> mpss(N);

  if(!input.restart)
    initialize_mpstates(mpos, mpss, fermionic::Quantum(input.N_elecs, input.N_spins), rgen, input.prefix, input.N_max_states);

  double esav = 1.0e8;

  cout << esav << "\t" << input.N_max_sweep_iter << endl;

  //
  // optimization with sweep algorithm
  //
  for(size_t iter = 0; iter < input.N_max_sweep_iter; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(mpos, mpss, input);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << endl;
    esav = eswp;
    if(std::fabs(edif) < input.tolerance) break;

  }

  std::ifstream fin("wave-mps_site-0.tmp");
  boost::archive::binary_iarchive iar(fin);
  iar >> mpss[0];

  for(int i = 1;i < N;++i){

     char name[50];

     sprintf(name,"right-mps_site-%d.tmp",i);

     std::ifstream fin2(name);
     boost::archive::binary_iarchive iar2(fin2);
     iar2 >> mpss[i];

  }

  mpsxx::save_mpx(mpss,"state");

  return esav;

}

