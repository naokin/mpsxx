#ifndef __MPSXX_DMRG_INPUT_H
#define __MPSXX_DMRG_INPUT_H

#include <vector>
#include <string>

namespace mpsxx {

enum DMRGALGORITHM { ONESITE, TWOSITE };

struct DMRGInput {

  bool
    restart;
  std::vector<double>
    energy;

  size_t
    N_sites;

  int
    N_spins;
  int
    N_elecs;

  size_t
    N_roots;
  int
    N_max_states;

  DMRGALGORITHM
    algorithm;

  size_t
    N_max_sweep_iter;
  double
    tolerance;

  std::string
    prefix;

  DMRGInput() {

    restart          = 0;

    N_sites          = 0;

    N_spins          = 0;
    N_elecs          = 0;

    N_roots          = 1;
    N_max_states     = 0;

    algorithm        = ONESITE;

    N_max_sweep_iter = 100;
    tolerance        = 1.0e-8;

    prefix           = ".";
  }
};

} // namespace mpsxx

#endif // __MPSXX_DMRG_INPUT_H
