#ifndef __MPSXX_COMPRESS_QC_MPOS_H
#define __MPSXX_COMPRESS_QC_MPOS_H

#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <MPX_types.h>
#include <symmetry/fermion.h>

namespace mpsxx {

size_t compress_mpos_2site (bool forward, MPO<double,fermion>& x, MPO<double,fermion>& y);

std::vector<size_t> compress_mpos_cycle (std::vector<MPO<double,fermion>>& mpos);

void compress_qc_mpos (
  const std::vector<int>& groups,
  const std::vector<std::vector<MPO<double,fermion>>>& mpos,
        std::vector<std::vector<MPO<double,fermion>>>& comp);

} // namespace mpsxx

#endif // __MPSXX_COMPRESS_QC_MPOS_H
