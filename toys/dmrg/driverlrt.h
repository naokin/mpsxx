#ifndef PROTOTYPE_DMRG_DRIVERLRT_H
#define PROTOTYPE_DMRG_DRIVERLRT_H

#include <vector>
#include <btas/DArray.h>

namespace lrt
{

void canonicalize              (bool forward,
                                const btas::DArray<4>& sdm1,
                                const btas::DArray<1>& val0,
                                const btas::DArray<3>& mps0,
                                const btas::DArray<3>& nul0,
                                      btas::DArray<3>& mps1);

void compute_sigma_vector      (const btas::DArray<4>& mpo0,
                                const btas::DArray<3>& mps0, const btas::DArray<3>& lstr0, const btas::DArray<3>& rstr0,
                                const btas::DArray<3>& mps1, const btas::DArray<3>& lstr1, const btas::DArray<3>& rstr1,
                                      btas::DArray<3>& sgv1);

std::vector<double> diagonalize(const btas::DArray<4>& mpo0,
                                const btas::DArray<3>& diag,
                                const btas::DArray<3>& lstr0,
                                const btas::DArray<3>& rstr0,
                                const btas::DArray<3>& mps0,
                                const std::vector< btas::DArray<3> >& lstr1,
                                const std::vector< btas::DArray<3> >& rstr1,
                                      std::vector< btas::DArray<3> >& mps1,
                                int nroot = 3, double tole = 1.0e-8, int max_iter = 100);

void compute_overlap_matrix(const std::vector< btas::DArray<3> >& wfns,
                            const std::vector< btas::DArray<3> >& lmps,
                            const std::vector< btas::DArray<3> >& rmps);

};

#endif // PROTOTYPE_DMRG_DRIVERLRT_H
