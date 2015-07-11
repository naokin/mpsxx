#ifndef _PROTOTYPE_DMRG_LRT_H
#define _PROTOTYPE_DMRG_LRT_H 1

#include <vector>

#include <legacy/DArray.h>
#include "mpsite.h"
#include "input.h"

namespace prototype
{

enum CALC_TYPE { TDA, RPA };

void RotateStorage
(const btas::DArray<2>& alpha, std::vector< btas::DArray<3> >& store, int nroots);

void RotateStorage
(const btas::DArray<2>& alphaRe, const btas::DArray<2>& alphaIm,
 std::vector< btas::DArray<3> >& storeRe, std::vector< btas::DArray<3> >& storeIm, int nroots);

void SolveCorrectionEquation
(bool forward, const MpSite& site0, MpSite& siteRe, MpSite& siteIm, const btas::DArray<3>& diag,
 const std::vector<double>& eigv, std::vector<double>& rnorm, int nroots, int mroots, int kroots, bool boundary, const CALC_TYPE& calc_type, double noise = 0.0);

void ComputeMatrixElements
(bool forward, const MpSite& site0, const MpSite& siteRe, const MpSite& siteIm,
 const std::vector<double>& eigv, btas::DArray<2>& a_subspace, btas::DArray<2>& b_subspace, btas::DArray<2>& s_subspace, btas::DArray<2>& d_subspace,
 std::vector<double>& bnorm, int nroots, int mroots, int kroots, bool boundary, bool last, const CALC_TYPE& calc_type);

void StorageInitialize(MpStorages& sites);

int lrt_sweep
(bool forward, MpStorages& sites0, MpStorages& sitesRe, MpStorages& sitesIm,
 btas::DArray<2>& a_subspace, btas::DArray<2>& b_subspace, btas::DArray<2>& s_subspace, btas::DArray<2>& d_subspace,
 std::vector<double>& eigv, std::vector<double>& rnorm, std::vector<double>& bnorm, int nroots, int mroots, int kroots, const CALC_TYPE& calc_type, double noise = 0.0);

std::vector<double> dmrg_lrt(const DmrgInput& input, MpStorages& sites);




void compute_derived_storage
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<3>& opr0,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket1,
                     btas::DArray<3>& fopr);

void transfer_derived_storage
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<3>& fopr,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket0,
                     btas::DArray<5>& gopr);

void transfer_derived_storage
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<5>& fopr,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket0,
                     btas::DArray<5>& gopr);

void compute_a_diagonal_elements
(bool forward, double e0,
               const btas::DArray<4>& mpo,
               const btas::DArray<3>& opr0,
               const btas::DArray<3>& copr,
               const btas::DArray<3>& bra1,
               const btas::DArray<3>& ket1,
                     btas::DArray<4>& aii);

void compute_a_nearest_elements
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<3>& fopr,
               const btas::DArray<3>& copr,
               const btas::DArray<3>& bra1,
               const btas::DArray<3>& ket0,
                     btas::DArray<4>& aij);

void compute_a_matrix_elements
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<5>& fopr,
               const btas::DArray<3>& copr,
               const btas::DArray<3>& bra1,
               const btas::DArray<3>& ket0,
                     btas::DArray<4>& aij);

void compute_b_nearest_elements
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<3>& fopr,
               const btas::DArray<3>& copr,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket1,
                     btas::DArray<4>& aij);

void compute_b_matrix_elements
(bool forward, const btas::DArray<4>& mpo,
               const btas::DArray<5>& fopr,
               const btas::DArray<3>& copr,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket1,
                     btas::DArray<4>& bij);

void compute_s_matrix_elements
(bool forward, const btas::DArray<2>& sopr,
               const btas::DArray<3>& bra1,
               const btas::DArray<3>& ket1,
                     btas::DArray<4>& sii);

void construct_effective_hamiltonian(MpStorages& sites, double eigenvalue);

void printing_symmetric_matrix(const btas::DArray<2>& a, int mcols = 8);

};

#endif // _PROTOTYPE_DMRG_LRT_H 1
