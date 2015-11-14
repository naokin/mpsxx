#ifndef __MPSXX_MAKE_SWEEP_H
#define __MPSXX_MAKE_SWEEP_H

#include <vector>

#include "input.h"

namespace mpsxx {

/// Perform 1 fwd./bwd. sweep iteration
double make_sweep(const std::vector<double>& E, const DMRGInput& input, const size_t& iroot = 0);

} // namespace mpsxx

#endif // __MPSXX_MAKE_SWEEP_H
