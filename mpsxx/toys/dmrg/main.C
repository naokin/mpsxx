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

  MpOprtrs mpos; /* hamiltonian in mpo rep.          */ /* L x (D, d, d, D)   */
  MpStates wfns; /* wavefunction for each site       */ /* L x (M, d, M)      */

  DnWeight lval; /* eigen values of left density     */ /* (L-1) x (4M)       */
  MpStates lmps; /* rotation matrix from left        */ /* (L-1) x (M, d, M)  */
  MpStates lnul; /* null-space of left density       */ /* (L-1) x (M, d, 3M) */
  Storages lstr; /* renormalized operator from left  */ /* L x (M, D, M)      */

  DnWeight rval; /* eigen values of right density    */ /* (L-1) x (4M)       */
  MpStates rmps; /* rotation matrix from right       */ /* (L-1) x (M, d, M)  */
  MpStates rnul; /* null-space of right density      */ /* (L-1) x (3M, d, M) */
  Storages rstr; /* renormalized operator from right */ /* L x (M, D, M)      */

  vector<MpStates> wavefuncs;
  vector<MpStates> lmpstates;
  vector<MpStates> rmpstates;

  //
  // set hamiltonian mpo
  //

  mpo_init(dmrginput, mpos);

  //
  // dmrg optimization
  //

  cout << "\t**************************************************"
       <<   "**************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::OPTIMIZATION"         << endl;
  cout << "\t**************************************************"
       <<   "**************************************************" << endl;
  dmrg   (cout, dmrginput, mpos, wfns, lval, lmps, lnul, lstr, rval, rmps, rnul, rstr);

  wavefuncs.push_back(wfns);
  lmpstates.push_back(lmps);
  rmpstates.push_back(rmps);

  cout << endl;

  cout << "\t**************************************************"
       <<   "**************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::OPTIMIZATION"         << endl;
  cout << "\t**************************************************"
       <<   "**************************************************" << endl;
  dmrg   (cout, dmrginput, mpos, wfns, lval, lmps, lnul, lstr, rval, rmps, rnul, rstr, wavefuncs, lmpstates, rmpstates);

  cout << endl;

  //
  // linear response for excited states
  //

//cout << "\t**************************************************"
//     <<   "**************************************************" << endl;
//cout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::LINEAR_RESPONSE"      << endl;
//cout << "\t**************************************************"
//     <<   "**************************************************" << endl;
//dmrglrt(cout, dmrginput, mpos, wfns, lval, lmps, lnul, lstr, rval, rmps, rnul, rstr, 4);

  cout << endl;

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
