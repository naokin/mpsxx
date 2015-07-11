#ifndef PROTOTYPE_DMRG_H
#define PROTOTYPE_DMRG_H

#include <iostream>
#include <vector>
#include <legacy/DArray.h>
#include "input.h"
#include "mpsite.h"

void mpo_init(const DmrgInput& input, MpStorages& sites);
void wfn_init(const DmrgInput& input, MpStorages& sites, int nroot = 0);
void str_init(const DmrgInput& input, MpStorages& sites, int nroot = 0);
double optimize(MpSite& site, const std::vector<double>& evals, double tole = 1.0e-8);
double dmrg(std::ostream& fout, const DmrgInput& input, MpStorages& sites, const std::vector<double>& evals);

#endif // PROTOTYPE_DMRG_H
