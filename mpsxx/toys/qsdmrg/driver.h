#ifndef _PROTOTYPE_DRIVER_H
#define _PROTOTYPE_DRIVER_H 1

#include <legacy/QSDArray.h>
#include <legacy/QSDblas.h>
#include <legacy/QSDlapack.h>

namespace prototype
{

void ComputeGuess
(bool forward, const btas::QSDArray<3>& mps0,
               const btas::QSDArray<3>& wfn0,
               const btas::QSDArray<3>& mps1,
                     btas::QSDArray<3>& wfn1);

void Canonicalize
(bool forward, const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& mps0, int M = 0);

void Canonicalize
(bool forward, const btas::QSDArray<4>& wfnx,
                     btas::QSDArray<3>& mps0,
                     btas::QSDArray<3>& wfn1, int M = 0);

void Renormalize
(bool forward, const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& opr0,
               const btas::QSDArray<3>& bra0,
               const btas::QSDArray<3>& ket0,
                     btas::QSDArray<3>& opr1);

void ComputeInverseGauge
(              const btas::QSDArray<2>& gauge_0,
                     btas::QSDArray<2>& gauge_i);

void Normalize
(                    btas::QSDArray<3>& wfn0);

void Orthogonalize
(              const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& wfn1);

void Orthogonalize
(bool forward, const btas::QSDArray<3>& mps0,
                     btas::QSDArray<3>& wfn1);


void ComputeDiagonal
(              const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
                     btas::QSDArray<3>& diag);

void ComputeDiagonal
(              const btas::QSDArray<4>& lmpo,
               const btas::QSDArray<4>& rmpo,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
                     btas::QSDArray<4>& diag);

void ComputeSigmaVector
(              const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
               const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& sgv0);

void ComputeSigmaVectorConj
(              const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
               const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& sgv0);

void ComputeSigmaVector
(              const btas::QSDArray<4>& lmpo,
               const btas::QSDArray<4>& rmpo,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
               const btas::QSDArray<4>& wfn0,
                     btas::QSDArray<4>& sgv0);

int ComputeOrthogonalTransform
(              const btas::DArray<2>& s_matrix,
                     btas::DArray<2>& u_matrix);

void TransformedMatrix
(              const btas::DArray<2>& h_matrix,
               const btas::DArray<2>& u_matrix,
                     btas::DArray<2>& h_transm);

int ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& s_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alpha);

int ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm);

void ComputeEigenvaluesTest
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace);

void AnalyzeTransferOperator
(              const btas::QSDArray<3>& bra,
               const btas::QSDArray<3>& ket);

};

#endif // _PROTOTYPE_DRIVER_H
