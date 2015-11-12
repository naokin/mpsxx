#ifndef PROTOTYPE_INPUT_H
#define PROTOTYPE_INPUT_H

#include <iostream>
#include <iomanip>

struct DmrgInput
{
  bool
    restart;

  int
    N_sites;
  int
    N_max_states;
  int
    N_phys_index;

  double
    J;
  double
    Jz;
  double
    hz;

  double
    tolerance;

  DmrgInput()
  : restart(0), N_sites(0), N_max_states(0), N_phys_index(2),
    J(1.0), Jz(1.0), hz(0.0), tolerance(1.0e-8)
  {
  }
};

std::ostream& operator<< (std::ostream& ost, const DmrgInput& input);

#endif // PROTOTYPE_INPUT_H
