#ifndef __MPSXX_DMRG_DENSE_OPTIMIZE_H
#define __MPSXX_DMRG_DENSE_OPTIMIZE_H

#include "MPX.h"

namespace mpsxx {

double optimize_1site (
        bool forward,
  const MPO  <double>& mpo0,
  const BLOCK<double>& lop0,
  const BLOCK<double>& rop0,
        MPS  <double>& mps0,
        MPS  <double>& mps1,
        BLOCK<double>& xop1);

double sweep (
  const MPOs<double>& mpos,
        MPSs<double>& mpss,
        std::vector<BLOCK<double>>& lops,
        std::vector<BLOCK<double>>& rops,
  const size_t& M);

double dmrg (
  const MPOs<double>& mpos,
        MPSs<double>& mpss,
  const size_t& M);

} // namespace mpsxx

#endif // __MPSXX_DMRG_DENSE_OPTIMIZE_H
