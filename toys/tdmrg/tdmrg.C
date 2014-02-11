#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <dense.h>
using namespace std;
using namespace btas;

#define DEBUG(msg) { cout << "DEBUG:" << msg << endl; }

int imag_tdmrg(ostream& fout, vector<DTensor<3> >& mpsts, bool restart,
               int L, int M, double J, double Jz, int N, double dt);
int real_tdmrg(ostream& fout, vector<ZTensor<3> >& mpsts, bool restart,
               int L, int M, double J, double Jz, double h, double f, int N, double dt);
int       dmrg(ostream& fout, vector<DTensor<3> >& mpsts, bool restart,
               int L, int M, double J, double Jz);

int main(int argc, char* argv[])
{
  // by default
  int    L  = 10;
  int    M  = 10;
  double J  = 1.0;
  double Jz = 1.0;
  double T  = 1.0;
  double dt = 0.01;
  double h  = 0.0;
  double f  = 0.0;

  string input;
  string output;
  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) input   = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) output  = argv[++iarg];
  }
  if(input.size() > 0) {
    ifstream fin(input.c_str());
    string entry;
    while(fin >> entry) {
      if(entry == "L" ) fin >> L;
      if(entry == "M" ) fin >> M;
      if(entry == "J" ) fin >> J;
      if(entry == "Jz") fin >> Jz;
      if(entry == "T" ) fin >> T;
      if(entry == "dt") fin >> dt;
      if(entry == "h" ) fin >> h;
      if(entry == "f" ) fin >> f;
    }
  }
  int N = (int) (T/dt + 0.5);

  streambuf *backup;
  backup = cout.rdbuf();

  ofstream fout;
  if(output.size() > 0) {
    fout.open(output.c_str());
    cout.rdbuf(fout.rdbuf());
  }

  vector<DTensor<3> > mpsts_g; // ground state mpss
  mpsts_g.reserve(L);
  // calling dmrg opt.
  cout << "DMRG SWEEP ALGO'S FOR SEARCHING GROUND STATE" << endl;
  dmrg(cout, mpsts_g, false, L, M, J, Jz);

//cout << "IMAGINARY TIME-EVOLUTION BY TEBD ALGO'S FOR GROUND STATE" << endl;
//imag_tdmrg(cout, mpsts_g, false, L, M, J, Jz, N, dt);
//imag_tdmrg(cout, mpsts_g, true, L, M, J, Jz, N, dt);

  vector<ZTensor<3> > mpsts_e; // time-evolved mpss
  mpsts_e.reserve(L);
  // copy real mps into complex mps
  for(int i = 0; i < mpsts_g.size(); ++i) {
    mpsts_e.push_back(ZTensor<3>(mpsts_g[i]));
  }
  // calling time-evolution
  cout << "REAL TIME-EVOLUTION BY TEBD ALGO'S FOR EXCITED STATE" << endl;
  real_tdmrg(cout, mpsts_e, true,  L, M, J, Jz, h, f, N, dt);
//real_tdmrg(cout, mpsts_e, false, L, M, J, Jz, h, f, N, dt);

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
