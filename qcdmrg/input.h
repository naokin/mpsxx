#ifndef _MPSXX_CXX11_DMRG_INPUT_H
#define _MPSXX_CXX11_DMRG_INPUT_H 1

#include <string>

namespace mpsxx {

enum DMRGALGORITHM { ONESITE, TWOSITE };

struct DmrgInput {

  bool
    restart;
  double
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

  DmrgInput() {

    restart          = 0;
    energy           = 0.0;

    N_sites          = 0;

    N_spins          = 0;
    N_elecs          = 0;

    N_roots          = 1;
    N_max_states     = 0;

    algorithm        = ONESITE;

    N_max_sweep_iter = 100;
    tolerance        = 1.0e-8;

    prefix           = "./";
  }
};

}; // namespace mpsxx

#endif // _MPSXX_CXX11_DMRG_INPUT_H
