#ifndef _PROTOTYPE_DMRG_H
#define _PROTOTYPE_DMRG_H 1

#include "mpsite.h"
#include "input.h"

namespace prototype
{

enum DMRG_ALGORITHM { ONESITE, TWOSITE };

namespace Heisenberg {

void construct_mpo(MpStorages& sites, const HeisenbergModel& params);

};

namespace Hubbard {

void construct_mpo(MpStorages& sites, const HubbardModel& params);

};

void initialize(MpStorages& sites, int M = 0);

double optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0, double noise = 0.0);

double optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0, double noise = 0.0);

double dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M = 0, double noise = 0.0);

double dmrg(const DmrgInput& input, MpStorages& sites);

};

#endif // _PROTOTYPE_DMRG_H
