#include <fstream>

#include "mpidefs.h"
#include "input.h"

void mpsxx::DMRGInput::parse (const std::string& fname)
{
  if(fname.size() == 0) return;

  Communicator world;

  if(world.rank() == 0) {
    std::ifstream fin(fname.c_str());
    std::string entry;
    while(fin >> entry) {
      if(entry == "restart")
        this->restart = true;
      if(entry == "N")
        fin >> this->N_sites;
      if(entry == "spin")
        fin >> this->N_spins;
      if(entry == "elec")
        fin >> this->N_elecs;
      if(entry == "M" || entry == "max_states")
        fin >> this->N_max_states;
      if(entry == "nroots")
        fin >> this->N_roots;
      if(entry == "tole" || entry == "tolerance")
        fin >> this->tolerance;
      if(entry == "noise")
        fin >> this->noise;
      if(entry == "onesite" || entry == "onedot")
        this->algorithm = mpsxx::ONESITE;
      if(entry == "twosite" || entry == "twodot")
        this->algorithm = mpsxx::TWOSITE;
      if(entry == "maxiter")
        fin >> this->N_max_sweep_iter;
    }
  }
#ifndef _SERIAL
  boost::mpi::broadcast(world,*this,0);
#endif
}
