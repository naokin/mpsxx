#ifndef __MPSXX_DMRG_DENSE_HUBBARD_H
#define __MPSXX_DMRG_DENSE_HUBBARD_H

#include <vector>

#include <legacy/DENSE/TArray.h>

#include "MPX.h"

namespace mpsxx {

/// construct MPOs of Hubbard model
void gen_hubbard_mpos (size_t N, MPOs<double>& mpos, double t, double U);

/// construct MPOs of Hubbard model having non-local single exponetial decaying terms
void gen_hubbard_mpos (size_t N, MPOs<double>& mpos, double t, double U, double F);

/// construct MPOs of Hubbard model having non-local single exponetial decaying terms
void gen_hubbard_mpos (size_t N, MPOs<double>& mpos, const btas::TArray<double,2>& t, const btas::TArray<double,2>& q, double U);

} // namespace mpsxx

#endif // __MPSXX_DMRG_DENSE_HUBBARD_H
