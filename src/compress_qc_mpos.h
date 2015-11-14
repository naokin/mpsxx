#ifndef __MPSXX_COMPRESS_QC_MPOS_H
#define __MPSXX_COMPRESS_QC_MPOS_H

#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <btas/QSPARSE/QSTArray.h>

#include "fermion.h"

namespace mpsxx {

size_t compress_mpos_2site (bool forward, btas::QSTArray<double,4,fermion>& x, btas::QSTArray<double,4,fermion>& y);

std::vector<size_t> compress_mpos_cycle (std::vector<btas::QSTArray<double,4,fermion>>& mpos);

/// Iterative MPOs merge
/// All pieces of MPOs must be in-core.
std::vector<int> compress_qc_mpos (
  const std::vector<int>& groups,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& mpos,
        std::vector<std::vector<btas::QSTArray<double,4,fermion>>>& comp);

/// Quantum# merge for single MPOs stored on disk.
void compress_qc_mpos (
  const size_t& N,
  const std::string& opname,
  const std::string& prefix = ".");

} // namespace mpsxx

#endif // __MPSXX_COMPRESS_QC_MPOS_H
