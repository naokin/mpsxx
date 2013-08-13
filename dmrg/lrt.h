#ifndef _PROTOTYPE_DMRG_LRT_H
#define _PROTOTYPE_DMRG_LRT_H 1

#include <vector>

#include <btas/DENSE/DArray.h>
#include <btas/QSPARSE/QSDArray.h>

#include "mpsite.h"
#include "input.h"

namespace prototype
{

enum CALC_TYPE { TDA, RPA };

void RotateStorage
(const btas::DArray<2>& alpha, std::vector< btas::QSDArray<3> >& store, int nroots);

void RotateStorage
(const btas::DArray<2>& alphaRe, const btas::DArray<2>& alphaIm,
 std::vector< btas::QSDArray<3> >& storeRe, std::vector< btas::QSDArray<3> >& storeIm, int nroots);

void SolveCorrectionEquation
(bool forward, const MpSite& site0, MpSite& siteRe, MpSite& siteIm, const btas::QSDArray<3>& diag,
 const std::vector<double>& eigv, std::vector<double>& rnorm, int nroots, int mroots, int kroots, bool boundary, const CALC_TYPE& calc_type);

void ComputeMatrixElements
(bool forward, const MpSite& site0, const MpSite& siteRe, const MpSite& siteIm,
 const std::vector<double>& eigv, btas::DArray<2>& a_subspace, btas::DArray<2>& b_subspace, btas::DArray<2>& s_subspace, btas::DArray<2>& d_subspace,
 std::vector<double>& bnorm, int nroots, int mroots, int kroots, bool boundary, bool last, const CALC_TYPE& calc_type);

void StorageInitialize(MpStorages& sites);

int lrt_sweep
(bool forward, MpStorages& sites0, MpStorages& sitesRe, MpStorages& sitesIm,
 btas::DArray<2>& a_subspace, btas::DArray<2>& b_subspace, btas::DArray<2>& s_subspace, btas::DArray<2>& d_subspace,
 std::vector<double>& eigv, std::vector<double>& rnorm, std::vector<double>& bnorm, int nroots, int mroots, int kroots, const CALC_TYPE& calc_type);

std::vector<double> dmrg_lrt(const DmrgInput& input, MpStorages& sites);

};

#endif // _PROTOTYPE_DMRG_H
