#ifndef __MPSXX_MPSGEN_MAKE_RANDOM_MPSS_H
#define __MPSXX_MPSGEN_MAKE_RANDOM_MPSS_H

#include <vector>
#include <btas/QSPARSE/QSTArray.h>

#include "input.h"
#include "fermion.h"

namespace mpsxx {

/// Generate quantum numbers for each boundary from MPO
///
/// \param sites quantum numbers for each site (physical index)
/// \param qt total quantum number of lattice
/// \param _max_quantum_blocks maximum number of quantum blocks on each boundary (0: no limitation by default)
std::vector<btas::Qshapes<fermion>>
generate_quantum_states (const std::vector<btas::Qshapes<fermion>>& sites, const fermion& qt, size_t max_quantum_blocks = 0);

void make_random_mpss (const DMRGInput& input, const size_t& iroot = 0);

}

#endif // __MPSXX_MPSGEN_MAKE_RANDOM_MPSS_H
