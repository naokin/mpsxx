#ifndef _PROTOTYPE_INPUT_H
#define _PROTOTYPE_INPUT_H 1

#include <iostream>
#include <iomanip>
#include <cstring>

namespace prototype
{

enum QUANTUM_MODELS
{
  HEISENBERG,
  HUBBARD
};

struct HeisenbergModel
{
  int    Nz;
  double J;
  double Jz;
  double Hz;

  HeisenbergModel(int _Nz = 1, double _J = 1.0, double _Jz = 1.0, double _Hz = 0.0)
  { Nz = _Nz; J = _J; Jz = _Jz; Hz = _Hz; }
};

struct HubbardModel
{
  double t;
  double U;

  HubbardModel(double _t = 1.0, double _U = 1.0)
  { t = _t; U = _U; }
};

struct DmrgInput
{
  bool
    restart;
  double
    energy;

  int
    N_spins;
  int
    N_elecs;

  int
    N_sites;
  int
    N_max_states;
  int
    N_roots;

  int
    N_max_sweep_iter;

  QUANTUM_MODELS
    model;

  HeisenbergModel
    heisenberg;

  HubbardModel
    hubbard;

  double
    tolerance;

  std::string
    prefix;

  DmrgInput() : restart(0),
                energy(0.0),
                N_spins(0),
                N_elecs(0),
                N_sites(0),
                N_max_states(0),
                N_roots(1),
                N_max_sweep_iter(100),
                model(HEISENBERG),
                tolerance(1.0e-8),
                prefix(".")
  { }
};

};

#endif // _PROTOTYPE_INPUT_H
