#include <iostream>
#include <vector>
#include <cstring>

#include "generate_hubbard_operators.h"

int main(int argc, char* argv[])
{
  using namespace mpsxx;

  bool enable_swap_sweep_dir = false;
  int iprint = 0;

  size_t N = 4;
  double t = -1.0;
  double u = 0.0;
  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-n") == 0) N = atoi(argv[++iarg]);
    if(strcmp(argv[iarg],"-t") == 0) t = atof(argv[++iarg]);
    if(strcmp(argv[iarg],"-u") == 0) u = atof(argv[++iarg]);
  }

  MPO<fermionic::Quantum> mpos(N);
  fermionic::generate_hubbard_operators(mpos, t, u);

  return 0;
}
