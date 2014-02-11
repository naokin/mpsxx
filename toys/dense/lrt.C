#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include <btas/DArray.h>
#include <btas/Dblas.h>
#include <btas/Dcontract.h>
#include "btas_template_specialize.h"

#include "lrt.h"
#include "driver.h"
#include "davidson.h"
using namespace btas;

extern double rgen();

void printing_matrix_elements(const DArray<2>& h_subspace, const DArray<2>& s_subspace)
{
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\t\t\tPRINTING SUBSPACE HAMILTONIAN" << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout.precision(8);
  for(int i = 0; i < h_subspace.rows(); ++i) {
    cout << "\t";
    for(int j = 0; j < h_subspace.cols(); ++j) {
      cout << setw(12) << fixed << h_subspace(i, j);
    }
    cout << endl;
  }
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\t\t\tPRINTING SUBSPACE OVERLAP" << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout.precision(8);
  for(int i = 0; i < s_subspace.rows(); ++i) {
    cout << "\t";
    for(int j = 0; j < s_subspace.cols(); ++j) {
      cout << setw(12) << fixed << s_subspace(i, j);
    }
    cout << endl;
  }
}

void printing_matrix_elements(const DArray<2>& a_subspace, const DArray<2>& b_subspace,
                              const DArray<2>& s_subspace, const DArray<2>& d_subspace)
{
  int nrows = a_subspace.rows();
  int ncols = a_subspace.cols();
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\t\t\tPRINTING SUBSPACE HAMILTONIAN" << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout.precision(8);
  for(int i = 0; i < nrows; ++i) {
    cout << "\t";
    for(int j = 0; j < ncols; ++j) {
      cout << setw(12) << fixed << a_subspace(i, j);
    }
    cout << " | ";
    for(int j = 0; j < ncols; ++j) {
      cout << setw(12) << fixed << b_subspace(i, j);
    }
    cout << endl;
  }
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\t\t\tPRINTING SUBSPACE OVERLAP" << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout.precision(8);
  for(int i = 0; i < nrows; ++i) {
    cout << "\t";
    for(int j = 0; j < ncols; ++j) {
      cout << setw(12) << fixed << s_subspace(i, j);
    }
    cout << " | ";
    for(int j = 0; j < ncols; ++j) {
      cout << setw(12) << fixed << d_subspace(i, j);
    }
    cout << endl;
  }
}

void checking_gauge_condition(const prototype::MpSite& lsite0, const prototype::MpSite& lsite1,
                              const prototype::MpSite& rsite0, const prototype::MpSite& rsite1, int nroots)
{
  cout.precision(16);
  for(int iroot = 1; iroot < nroots; ++iroot) {
    DArray<4> LHS; // LHS = C(I)R(0) + C(0)R(I)
    Dgemm(NoTrans, NoTrans, 1.0, lsite1.wfnc[iroot], rsite0.rmps[0], 1.0, LHS);
    Dgemm(NoTrans, NoTrans, 1.0, lsite0.wfnc[0], rsite1.rmps[iroot], 1.0, LHS);

    DArray<4> RHS; // RHS = L(I)C(0) + L(0)C(I)
    Dgemm(NoTrans, NoTrans, 1.0, lsite1.lmps[iroot], rsite0.wfnc[0], 1.0, RHS);
    Dgemm(NoTrans, NoTrans, 1.0, lsite0.lmps[0], rsite1.wfnc[iroot], 1.0, RHS);

    Daxpy(-1.0, RHS, LHS); // LHS = LHS - RHS
    cout << "\tChecking gauge-condition[" << setw(2) << iroot << "]: subnorm = " << setw(20) << fixed << Ddot(LHS, LHS) << endl;
  }
}

void prototype::RotateStorage
(const DArray<2>& alpha, vector< DArray<3> >& store, int nroots)
{
  int mroots = 1 + alpha.rows();
  assert(store.size() >= mroots);
  vector< DArray<3> > tmpstr(nroots);
  for(int iroot = 1; iroot < nroots; ++iroot)
    for(int jroot = 1; jroot < mroots; ++jroot)
      Daxpy(alpha(jroot-1, iroot-1), store[jroot], tmpstr[iroot]);

  for(int iroot = 1; iroot < nroots; ++iroot)
    Dcopy(tmpstr[iroot], store[iroot]);
}

void prototype::RotateStorage
(const DArray<2>& alphaRe, const DArray<2>& alphaIm, vector< DArray<3> >& storeRe, vector< DArray<3> >& storeIm, int nroots)
{
  int mroots = 1 + alphaRe.rows();
  assert(storeRe.size() >= mroots);
  assert(storeIm.size() >= mroots);
  vector< DArray<3> > tmpstrRe(nroots);
  vector< DArray<3> > tmpstrIm(nroots);
  for(int iroot = 1; iroot < nroots; ++iroot) {
    for(int jroot = 1; jroot < mroots; ++jroot) {
      Daxpy(alphaRe(jroot-1, iroot-1), storeRe[jroot], tmpstrRe[iroot]);
      Daxpy(alphaIm(jroot-1, iroot-1), storeIm[jroot], tmpstrRe[iroot]);
      Daxpy(alphaRe(jroot-1, iroot-1), storeIm[jroot], tmpstrIm[iroot]);
      Daxpy(alphaIm(jroot-1, iroot-1), storeRe[jroot], tmpstrIm[iroot]);
    }
  }

  for(int iroot = 1; iroot < nroots; ++iroot) {
    Dcopy(tmpstrRe[iroot], storeRe[iroot]);
    Dcopy(tmpstrIm[iroot], storeIm[iroot]);
  }
}

void prototype::SolveCorrectionEquation
(bool forward, const MpSite& site0, MpSite& siteRe, MpSite& siteIm, const DArray<3>& diag,
 const vector<double>& eigv, vector<double>& rnorm, int nroots, int mroots, int kroots, bool boundary, const CALC_TYPE& calc_type, double noise)
{
  if(mroots == 1) {
    DArray<3> tmpwfc(site0.wfnc[0]);
    for(int iroot = 1; iroot < nroots; ++iroot) {
      siteRe.wfnc[iroot].free();
      ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], tmpwfc, siteRe.wfnc[iroot]);

      Orthogonalize(site0.wfnc[0], siteRe.wfnc[iroot]);
      for(int jroot = 1; jroot < iroot; ++jroot) {
        Normalize(siteRe.wfnc[iroot]);
        Orthogonalize(site0.wfnc[0], siteRe.wfnc[iroot]);
      }
      Normalize(siteRe.wfnc[iroot]);
      Dcopy(siteRe.wfnc[iroot], tmpwfc);

      Dcopy(site0.wfnc[0], siteIm.wfnc[iroot]);
      siteIm.wfnc[iroot] = 0.0;
// DEBUG FOR TEST
//    siteRe.wfnc[iroot] = rgen;
//    Normalize(siteRe.wfnc[iroot]);
//    Orthogonalize(site0.wfnc[0], siteRe.wfnc[iroot]);

//    siteIm.wfnc[iroot] = rgen;
//    Normalize(siteIm.wfnc[iroot]);
//    Orthogonalize(site0.wfnc[0], siteIm.wfnc[iroot]);
    }

    return;
  }

  vector< DArray<3> > sgvxRe(nroots+kroots-1);
  vector< DArray<3> > sgvxIm(nroots+kroots-1);

  // diagonal block
  for(int iroot = 1; iroot < nroots+kroots-1; ++iroot) {
    ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], siteRe.wfnc[iroot], sgvxRe[iroot]);
    Daxpy(-eigv[0], siteRe.wfnc[iroot], sgvxRe[iroot]);
    ComputeSigmaVector(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxRe[iroot]);
    ComputeSigmaVector(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgvxRe[iroot]);
  }
  for(int iroot = 1; iroot < nroots+kroots-1; ++iroot) {
    ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], siteIm.wfnc[iroot], sgvxIm[iroot]);
    Daxpy(-eigv[0], siteIm.wfnc[iroot], sgvxIm[iroot]);
    ComputeSigmaVector(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxIm[iroot]);
    ComputeSigmaVector(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgvxIm[iroot]);
  }
  // off-diagonal block
  if(calc_type == RPA) {

  for(int iroot = 1; iroot < nroots+kroots-1; ++iroot) {
    ComputeSigmaVectorConj(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxRe[iroot]);
    ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgvxRe[iroot]);
  }
  for(int iroot = 1; iroot < nroots+kroots-1; ++iroot) {
    ComputeSigmaVectorConj(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxIm[iroot]);
    ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgvxIm[iroot]);
  }

  }

  // residual norm & correction
  for(int iroot = 1; iroot < nroots+kroots-1; ++iroot) {
    DArray<3> errvRe(sgvxRe[iroot]);
    Daxpy(-eigv[iroot], siteRe.wfnc[iroot], errvRe);

    DArray<3> errvIm(sgvxIm[iroot]);
    Daxpy(+eigv[iroot], siteIm.wfnc[iroot], errvIm);

    double normRe = Ddot(errvRe, errvRe);
    double normIm = Ddot(errvIm, errvIm);

    if(iroot < nroots)
      rnorm[iroot] += (normRe + normIm);

    if(iroot < kroots) continue;

    davidson::precondition(+eigv[iroot], diag, errvRe);
//  util::Normalize(errvRe);
    util::Randomize(noise, errvRe, rgen);
    util::Normalize(errvRe);
    Orthogonalize(site0.wfnc[0], errvRe);

    if(calc_type == RPA) {
      davidson::precondition(-eigv[iroot], diag, errvIm);
//    util::Normalize(errvIm);
      util::Randomize(noise, errvIm, rgen);
      util::Normalize(errvIm);
      Orthogonalize(site0.wfnc[0], errvIm);
    }
    else {
      errvIm = 0.0;
    }

    Dcopy(errvRe, siteRe.wfnc[mroots+iroot-kroots]);
    Dcopy(errvIm, siteIm.wfnc[mroots+iroot-kroots]);
  }
}

void prototype::ComputeMatrixElements
(bool forward, const MpSite& site0, const MpSite& siteRe, const MpSite& siteIm,
 const vector<double>& eigv, DArray<2>& a_subspace, DArray<2>& b_subspace, DArray<2>& s_subspace, DArray<2>& d_subspace,
       vector<double>& bnorm, int nroots, int mroots, int kroots, bool boundary, bool last, const CALC_TYPE& calc_type)
{
  int lroots = mroots + nroots - 1;

  // projection to orthogonal subspace
  vector< DArray<3> > qfncRe(lroots);
  vector< DArray<3> > qfncIm(lroots);
  for(int iroot = 1; iroot < lroots; ++iroot) {
    Dcopy(siteRe.wfnc[iroot], qfncRe[iroot]);
    Dcopy(siteIm.wfnc[iroot], qfncIm[iroot]);

    if(last) continue;

    if(forward) {
      Orthogonalize(1, site0.lmps[0], qfncRe[iroot]);
      Orthogonalize(1, site0.lmps[0], qfncIm[iroot]);
    }
    else {
      Orthogonalize(0, site0.rmps[0], qfncRe[iroot]);
      Orthogonalize(0, site0.rmps[0], qfncIm[iroot]);
    }
  }

  // Xi Sij Xj
  for(int iroot = 1; iroot < lroots; ++iroot) {
    double sii = Ddot(qfncRe[iroot], qfncRe[iroot]);
    s_subspace(iroot-1, iroot-1) += sii;
    for(int jroot = 1; jroot < iroot; ++jroot) {
      double sij = Ddot(qfncRe[iroot], qfncRe[jroot]);
      s_subspace(iroot-1, jroot-1) += sij;
      s_subspace(jroot-1, iroot-1) += sij;
    }
  }
  // Yi Sij Yj
  for(int iroot = 1; iroot < lroots; ++iroot) {
    double sii = Ddot(qfncIm[iroot], qfncIm[iroot]);
    bnorm[iroot] += sii; // save Y-norm to compute RPA energy
    s_subspace(iroot-1, iroot-1) -= sii;
    for(int jroot = 1; jroot < iroot; ++jroot) {
      double sij = Ddot(qfncIm[iroot], qfncIm[jroot]);
      s_subspace(iroot-1, jroot-1) -= sij;
      s_subspace(jroot-1, iroot-1) -= sij;
    }
  }
  // Xi Sij Yj
  for(int iroot = 1; iroot < lroots; ++iroot) {
    for(int jroot = 1; jroot < iroot; ++jroot) {
      double dij = Ddot(qfncRe[iroot], qfncIm[jroot]);
      double dji = Ddot(qfncRe[jroot], qfncIm[iroot]);
      d_subspace(iroot-1, jroot-1) += (dij - dji);
      d_subspace(jroot-1, iroot-1) += (dji - dij);
    }
  }
  // compute reduced A-matrix & B-matrix / 1-st order terms
  // Aij += Xi Hij Xj + Yi Hij Yj, Bij = Xi Hij Yj + Yi Hij Xj
  {
    vector< DArray<3> > sgvxRe(lroots);
    vector< DArray<3> > sgv0Re(lroots);

    vector< DArray<3> > sgvxIm(lroots);
    vector< DArray<3> > sgv0Im(lroots);

    for(int iroot = 1; iroot < lroots; ++iroot) {
        ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], qfncRe[iroot], sgvxRe[iroot]);
        Daxpy(-eigv[0], qfncRe[iroot], sgvxRe[iroot]);
      if(forward)
        ComputeSigmaVector(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Re[iroot]);
      else
        ComputeSigmaVector(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgv0Re[iroot]);
    }
    for(int iroot = 1; iroot < lroots; ++iroot) {
        ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], qfncIm[iroot], sgvxIm[iroot]);
        Daxpy(-eigv[0], qfncIm[iroot], sgvxIm[iroot]);
      if(forward)
        ComputeSigmaVector(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Im[iroot]);
      else
        ComputeSigmaVector(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgv0Im[iroot]);
    }

//cout.precision(8);
//cout << "debug: check a_subspace (real x real)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = Ddot(qfncRe[iroot], sgvxRe[iroot])
                 + Ddot(qfncRe[iroot], sgv0Re[iroot]) * 2.0;
      a_subspace(iroot-1, iroot-1) += vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = Ddot(qfncRe[iroot], sgvxRe[jroot])
                   + Ddot(qfncRe[iroot], sgv0Re[jroot])
                   + Ddot(qfncRe[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        a_subspace(iroot-1, jroot-1) += vij;
        a_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << vii << endl;
    }
//cout << "debug: check a_subspace (imag x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = Ddot(qfncIm[iroot], sgvxIm[iroot])
                 + Ddot(qfncIm[iroot], sgv0Im[iroot]) * 2.0;
      a_subspace(iroot-1, iroot-1) += vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = Ddot(qfncIm[iroot], sgvxIm[jroot])
                   + Ddot(qfncIm[iroot], sgv0Im[jroot])
                   + Ddot(qfncIm[jroot], sgv0Im[iroot]);
//cout << setw(12) << vij;
        a_subspace(iroot-1, jroot-1) += vij;
        a_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << vii << endl;
    }
//cout << "debug: check b_subspace (real x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = Ddot(qfncRe[iroot], sgvxIm[iroot])
                 + Ddot(qfncRe[iroot], sgv0Im[iroot])
                 + Ddot(qfncIm[iroot], sgv0Re[iroot]);
      b_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = Ddot(qfncRe[iroot], sgvxIm[jroot])
                   + Ddot(qfncRe[iroot], sgv0Im[jroot])
                   + Ddot(qfncRe[jroot], sgv0Im[iroot])
                   + Ddot(qfncIm[iroot], sgvxRe[jroot])
                   + Ddot(qfncIm[iroot], sgv0Re[jroot])
                   + Ddot(qfncIm[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        b_subspace(iroot-1, jroot-1) += vij;
        b_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << 2.0 * vii << endl;
    }
  }

//cout << "\t****************************************************************************************************" << endl; 
//cout << "\tdebug @ ComputeMatrixElements: after computed diagonal-block" << endl;
//printing_matrix_elements(a_subspace, b_subspace, s_subspace, d_subspace);
//cout << "\t****************************************************************************************************" << endl; 

// DEBUG FOR TEST
//cout << "debug: suppose B = 0" << endl;
//return;

  if(boundary || calc_type != RPA) return;

  DArray<2> gauge0;
  if(forward)
    Dgemm(NoTrans, ConjTrans, 1.0, site0.wfnc[0], site0.rmps[0], 1.0, gauge0);
  else
    Dgemm(ConjTrans, NoTrans, 1.0, site0.lmps[0], site0.wfnc[0], 1.0, gauge0);

  // compute reduced A-matrix & B-matrix / 2-nd order terms
  // Aij += Xi Wij Yj + Yi Wij Xj, Bij = Xi Wij Xj + Yi Wij Yj
  {
    vector< DArray<3> > sgv0Re(lroots);
    vector< DArray<3> > sgv0Im(lroots);

    for(int iroot = 1; iroot < lroots; ++iroot) {
      qfncRe[iroot].free();
      if(forward) {
        Dgemm(NoTrans, NoTrans, 1.0, gauge0, siteRe.rmps[iroot], 1.0, qfncRe[iroot]);
        Orthogonalize(site0.wfnc[0], qfncRe[iroot]);
        ComputeSigmaVectorConj(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Re[iroot]);
//cout << "debug: forward: qfncRe[" << iroot << "] = " << qfncRe[iroot] << endl;
      }
      else {
        Dgemm(NoTrans, NoTrans, 1.0, siteRe.lmps[iroot], gauge0, 1.0, qfncRe[iroot]);
        Orthogonalize(site0.wfnc[0], qfncRe[iroot]);
        ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgv0Re[iroot]);
//cout << "debug: backward: qfncRe[" << iroot << "] = " << qfncRe[iroot] << endl;
      }
    }
    for(int iroot = 1; iroot < lroots; ++iroot) {
      qfncIm[iroot].free();
      if(forward) {
        Dgemm(NoTrans, NoTrans, 1.0, gauge0, siteIm.rmps[iroot], 1.0, qfncIm[iroot]);
        Orthogonalize(site0.wfnc[0], qfncIm[iroot]);
        ComputeSigmaVectorConj(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Im[iroot]);
//cout << "debug: forward: qfncIm[" << iroot << "] = " << qfncIm[iroot] << endl;
      }
      else {
        Dgemm(NoTrans, NoTrans, 1.0, siteIm.lmps[iroot], gauge0, 1.0, qfncIm[iroot]);
        Orthogonalize(site0.wfnc[0], qfncIm[iroot]);
        ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgv0Im[iroot]);
//cout << "debug: backward: qfncIm[" << iroot << "] = " << qfncIm[iroot] << endl;
      }
    }

//cout << "debug: check a_subspace (real x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = Ddot(qfncRe[iroot], sgv0Im[iroot])
                 + Ddot(qfncIm[iroot], sgv0Re[iroot]);
      a_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = Ddot(qfncRe[iroot], sgv0Im[jroot])
                   + Ddot(qfncRe[jroot], sgv0Im[iroot])
                   + Ddot(qfncIm[iroot], sgv0Re[jroot])
                   + Ddot(qfncIm[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        a_subspace(iroot-1, jroot-1) += vij;
        a_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << 2.0 * vii << endl;
    }
//cout << "debug: check b_subspace (real x real)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = Ddot(qfncRe[iroot], sgv0Re[iroot]);
      b_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = Ddot(qfncRe[iroot], sgv0Re[jroot])
                   + Ddot(qfncRe[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        b_subspace(iroot-1, jroot-1) += vij;
        b_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << 2.0 * vii << endl;
    }
//cout << "debug: check b_subspace (imag x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = Ddot(qfncIm[iroot], sgv0Im[iroot]);
      b_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = Ddot(qfncIm[iroot], sgv0Im[jroot])
                   + Ddot(qfncIm[jroot], sgv0Im[iroot]);
//cout << setw(12) << vij;
        b_subspace(iroot-1, jroot-1) += vij;
        b_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << 2.0 * vii << endl;
    }
  }
}

void prototype::StorageInitialize(MpStorages& sites)
{
  int L = sites.size();

  for(int i = 0; i < L; ++i) {
    for(int iroot = 0; iroot < sites[i].nstore; ++iroot) {
      sites[i].lopr[iroot].resize(1, 1, 1);
      sites[i].lopr[iroot] = 0.0;
      sites[i].ropr[iroot].resize(1, 1, 1);
      sites[i].ropr[iroot] = 0.0;
    }
    sites[i].save(i);
  }
}

int prototype::lrt_sweep(bool forward, MpStorages& sites0, MpStorages& sitesRe, MpStorages& sitesIm,
                         DArray<2>& a_subspace, DArray<2>& b_subspace, DArray<2>& s_subspace, DArray<2>& d_subspace,
                         vector<double>& eigv, vector<double>& rnorm, vector<double>& bnorm, int nroots, int mroots, int kroots, const CALC_TYPE& calc_type, double noise)
{
  DArray<2> alphaRe;
  DArray<2> alphaIm;

  if(mroots > 1) {
    printing_matrix_elements(a_subspace, b_subspace, s_subspace, d_subspace);

    DArray<1> eigvtmp;
    int reduced_mroots = 1 + ComputeEigenvalues(a_subspace, b_subspace, s_subspace, d_subspace, eigvtmp, alphaRe, alphaIm);
// FOR TEST
//  TEST::ComputeEigenvalues(a_subspace, b_subspace, s_subspace, d_subspace, eigvtmp, alphaRe, alphaIm);
// END TEST

    if(reduced_mroots < mroots)
      mroots = reduced_mroots;

    cout << "\tEigenvalues" << endl;
    for(int i = 1; i < nroots+kroots-1; ++i) {
      eigv[i] = eigvtmp(i-1);
      cout << "\t\tState [ " << setw(2) << i << " ]: " << setw(16) << fixed << setprecision(12) << eigv[i] << endl;
    }
  }

  int lroots = mroots + nroots - 1;

  int L = sites0.size();

  a_subspace.resize(lroots-1, lroots-1); a_subspace = 0.0;
  b_subspace.resize(lroots-1, lroots-1); b_subspace = 0.0;
  s_subspace.resize(lroots-1, lroots-1); s_subspace = 0.0;
  d_subspace.resize(lroots-1, lroots-1); d_subspace = 0.0;

  rnorm.resize(nroots); fill(rnorm.begin(), rnorm.end(), 0.0);
  bnorm.resize(lroots); fill(bnorm.begin(), bnorm.end(), 0.0);

//if(forward)
  {

  // fowrad sweep

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\t" << "FORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\tnroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << ", lroots = " << lroots << endl;
  cout << "\t====================================================================================================" << endl;

  sites0 [0].load(0);
  sitesRe[0].load(0);
  sitesIm[0].load(0);

  RotateStorage(alphaRe, alphaIm, sitesRe[0].wfnc, sitesIm[0].wfnc, mroots);

  bool boundary = true;

  for(int i = 0; i < L-1; ++i) {
    cout << "\t\t\tSTARTING BLOCK ITERATION FOR SITE " << i << endl;
    cout << "\t====================================================================================================" << endl;
    sites0 [i+1].load(i+1);
    sitesRe[i+1].load(i+1);
    sitesIm[i+1].load(i+1);

    DArray<3> diag;
    ComputeDiagonal(sites0[i].mpo, sites0[i].lopr[0], sites0[i].ropr[0], diag);
    for(DArray<3>::iterator it = diag.begin(); it != diag.end(); ++it) *it -= eigv[0];

    RotateStorage(alphaRe, alphaIm, sitesRe[i].ropr, sitesIm[i].ropr, mroots);

    cout << "\tSolving correction equation @ site " << i << endl;
    SolveCorrectionEquation(1, sites0[i], sitesRe[i], sitesIm[i], diag, eigv, rnorm, nroots, mroots, kroots, boundary, calc_type, noise);

    RotateStorage(alphaRe, alphaIm, sitesRe[i+1].rmps, sitesIm[i+1].rmps, mroots);

    for(int iroot = 1; iroot < mroots; ++iroot) {
      sitesRe[i+1].wfnc[iroot].free();
      ComputeGuess(1, sites0[i].lmps[0], sites0[i].wfnc[0], sitesRe[i+1].rmps[iroot], sitesRe[i+1].wfnc[iroot]);
      ComputeGuess(1, sites0[i].lmps[0], sitesRe[i].wfnc[iroot], sites0[i+1].rmps[0], sitesRe[i+1].wfnc[iroot]);
      Orthogonalize(sites0[i+1].wfnc[0], sitesRe[i+1].wfnc[iroot]);

      sitesIm[i+1].wfnc[iroot].free();
      ComputeGuess(1, sites0[i].lmps[0], sites0[i].wfnc[0], sitesIm[i+1].rmps[iroot], sitesIm[i+1].wfnc[iroot]);
      ComputeGuess(1, sites0[i].lmps[0], sitesIm[i].wfnc[iroot], sites0[i+1].rmps[0], sitesIm[i+1].wfnc[iroot]);
      Orthogonalize(sites0[i+1].wfnc[0], sitesIm[i+1].wfnc[iroot]);
    }
    for(int iroot = mroots; iroot < lroots; ++iroot) {
      sitesRe[i+1].wfnc[iroot].free();
      ComputeGuess(1, sites0[i].lmps[0], sitesRe[i].wfnc[iroot], sites0[i+1].rmps[0], sitesRe[i+1].wfnc[iroot]);
      Orthogonalize(sites0[i+1].wfnc[0], sitesRe[i+1].wfnc[iroot]);

      sitesIm[i+1].wfnc[iroot].free();
      ComputeGuess(1, sites0[i].lmps[0], sitesIm[i].wfnc[iroot], sites0[i+1].rmps[0], sitesIm[i+1].wfnc[iroot]);
      Orthogonalize(sites0[i+1].wfnc[0], sitesIm[i+1].wfnc[iroot]);
    }

    boundary = false;

    DArray<2> gauge_0;
    Dgemm(ConjTrans, NoTrans, 1.0, sites0[i].lmps[0], sites0[i].wfnc[0], 1.0, gauge_0);
    DArray<2> gauge_i;
    ComputeInverseGauge(gauge_0, gauge_i);

//DArray<2> check0i;
//Dgemm(NoTrans, NoTrans, 1.0, gauge_0, gauge_i, 1.0, check0i);
//cout << "debug: gauge-trans 0-i check = " << check0i << endl;
//DArray<2> checki0;
//Dgemm(NoTrans, NoTrans, 1.0, gauge_i, gauge_0, 1.0, checki0);
//cout << "debug: gauge-trans i-0 check = " << checki0 << endl;

    for(int iroot = 1; iroot < lroots; ++iroot) {
      DArray<3> tmpRe(sitesRe[i].wfnc[iroot]);
      Orthogonalize(1, sites0[i].lmps[0], tmpRe);
      sitesRe[i].lmps[iroot].free();
      Dgemm(NoTrans, NoTrans, 1.0, tmpRe, gauge_i, 1.0, sitesRe[i].lmps[iroot]);

      sitesRe[i+1].lopr[iroot].free();
      Renormalize(1, sites0[i].mpo, sitesRe[i].lopr[iroot], sites0[i].lmps[0], sites0[i].lmps[0], sitesRe[i+1].lopr[iroot]);
      Renormalize(1, sites0[i].mpo, sites0[i].lopr[0], sites0[i].lmps[0], sitesRe[i].lmps[iroot], sitesRe[i+1].lopr[iroot]);

      DArray<3> tmpIm(sitesIm[i].wfnc[iroot]);
      Orthogonalize(1, sites0[i].lmps[0], tmpIm);
      sitesIm[i].lmps[iroot].free();
      Dgemm(NoTrans, NoTrans, 1.0, tmpIm, gauge_i, 1.0, sitesIm[i].lmps[iroot]);

      sitesIm[i+1].lopr[iroot].free();
      Renormalize(1, sites0[i].mpo, sitesIm[i].lopr[iroot], sites0[i].lmps[0], sites0[i].lmps[0], sitesIm[i+1].lopr[iroot]);
      Renormalize(1, sites0[i].mpo, sites0[i].lopr[0], sites0[i].lmps[0], sitesIm[i].lmps[iroot], sitesIm[i+1].lopr[iroot]);
    }

    sites0 [i].save(i);
    sitesRe[i].save(i);
    sitesIm[i].save(i);
    cout << "\t====================================================================================================" << endl;
  }

  {
    DArray<3> diag;
    ComputeDiagonal(sites0[L-1].mpo, sites0[L-1].lopr[0], sites0[L-1].ropr[0], diag);
    for(DArray<3>::iterator it = diag.begin(); it != diag.end(); ++it) *it -= eigv[0];

    RotateStorage(alphaRe, alphaIm, sitesRe[L-1].ropr, sitesIm[L-1].ropr, mroots);

    cout << "\tSolving correction equation @ site " << L-1 << endl;
    SolveCorrectionEquation(1, sites0[L-1], sitesRe[L-1], sitesIm[L-1], diag, eigv, rnorm, nroots, mroots, kroots, boundary, calc_type, noise);

    sites0 [L-1].save(L-1);
    sitesRe[L-1].save(L-1);
    sitesIm[L-1].save(L-1);
  }

  } // if(forward)
//else
  {

  // backward sweep

  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\t" << "BACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\tnroots = " << nroots << ", mroots = " << mroots << ", kroots = " << kroots << ", lroots = " << lroots << endl;
  cout << "\t====================================================================================================" << endl;

  sites0 [L-1].load(L-1);
  sitesRe[L-1].load(L-1);
  sitesIm[L-1].load(L-1);

  bool boundary = true;

  for(int i = L-1; i > 0; --i) {
    cout << "\t\t\tSTARTING BLOCK ITERATION FOR SITE " << i << endl;
    cout << "\t====================================================================================================" << endl;
    sites0 [i-1].load(i-1);
    sitesRe[i-1].load(i-1);
    sitesIm[i-1].load(i-1);

    for(int iroot = 1; iroot < lroots; ++iroot) {
      sitesRe[i-1].wfnc[iroot].free();
      ComputeGuess(0, sites0[i].rmps[0], sites0[i].wfnc[0], sitesRe[i-1].lmps[iroot], sitesRe[i-1].wfnc[iroot]);
      ComputeGuess(0, sites0[i].rmps[0], sitesRe[i].wfnc[iroot], sites0[i-1].lmps[0], sitesRe[i-1].wfnc[iroot]);
      Orthogonalize(sites0[i-1].wfnc[0], sitesRe[i-1].wfnc[iroot]);

      sitesIm[i-1].wfnc[iroot].free();
      ComputeGuess(0, sites0[i].rmps[0], sites0[i].wfnc[0], sitesIm[i-1].lmps[iroot], sitesIm[i-1].wfnc[iroot]);
      ComputeGuess(0, sites0[i].rmps[0], sitesIm[i].wfnc[iroot], sites0[i-1].lmps[0], sitesIm[i-1].wfnc[iroot]);
      Orthogonalize(sites0[i-1].wfnc[0], sitesIm[i-1].wfnc[iroot]);
    }

    cout << "\tComputing matrix elements @ site " << i << endl;
    ComputeMatrixElements(0, sites0[i], sitesRe[i], sitesIm[i], eigv, a_subspace, b_subspace, s_subspace, d_subspace, bnorm, nroots, mroots, kroots, boundary, false, calc_type);
//  printing_matrix_elements(a_subspace, b_subspace, s_subspace, d_subspace);

    boundary = false;

    DArray<2> gauge_0;
    Dgemm(NoTrans, ConjTrans, 1.0, sites0[i].wfnc[0], sites0[i].rmps[0], 1.0, gauge_0);
    DArray<2> gauge_i;
    ComputeInverseGauge(gauge_0, gauge_i);

    for(int iroot = 1; iroot < lroots; ++iroot) {
      DArray<3> tmpRe(sitesRe[i].wfnc[iroot]);
      Orthogonalize(0, sites0[i].rmps[0], tmpRe);
      sitesRe[i].rmps[iroot].free();
      Dgemm(NoTrans, NoTrans, 1.0, gauge_i, tmpRe, 1.0, sitesRe[i].rmps[iroot]);

      sitesRe[i-1].ropr[iroot].free();
      Renormalize(0, sites0[i].mpo, sitesRe[i].ropr[iroot], sites0[i].rmps[0], sites0[i].rmps[0], sitesRe[i-1].ropr[iroot]);
      Renormalize(0, sites0[i].mpo, sites0[i].ropr[0], sites0[i].rmps[0], sitesRe[i].rmps[iroot], sitesRe[i-1].ropr[iroot]);

      DArray<3> tmpIm(sitesIm[i].wfnc[iroot]);
      Orthogonalize(0, sites0[i].rmps[0], tmpIm);
      sitesIm[i].rmps[iroot].free();
      Dgemm(NoTrans, NoTrans, 1.0, gauge_i, tmpIm, 1.0, sitesIm[i].rmps[iroot]);

      sitesIm[i-1].ropr[iroot].free();
      Renormalize(0, sites0[i].mpo, sitesIm[i].ropr[iroot], sites0[i].rmps[0], sites0[i].rmps[0], sitesIm[i-1].ropr[iroot]);
      Renormalize(0, sites0[i].mpo, sites0[i].ropr[0], sites0[i].rmps[0], sitesIm[i].rmps[iroot], sitesIm[i-1].ropr[iroot]);
    }

    sites0 [i].save(i);
    sitesRe[i].save(i);
    sitesIm[i].save(i);
    cout << "\t====================================================================================================" << endl;
  }

  {
    cout << "\tComputing matrix elements @ site " << 0 << endl;
    ComputeMatrixElements(0, sites0[0], sitesRe[0], sitesIm[0], eigv, a_subspace, b_subspace, s_subspace, d_subspace, bnorm, nroots, mroots, kroots, true, true, calc_type);
//  printing_matrix_elements(a_subspace, b_subspace, s_subspace, d_subspace);

    sites0 [0].save(0);
    sitesRe[0].save(0);
    sitesIm[0].save(0);
  }

  } // if(forward) else

  return lroots;
}

vector<double> prototype::dmrg_lrt(const DmrgInput& input, MpStorages& sites0)
{
  const int L = input.N_sites;
  const int M = input.N_max_states;

  const int nroots = input.N_roots;

  int mxritz = max(2+4*(nroots-1), 20);
  int nstore = (nroots == 1) ? 1 : (nroots+mxritz);

  MpStorages sitesRe(L, MpSite("state-Re", nstore, input.prefix));
  StorageInitialize(sitesRe);

  MpStorages sitesIm(L, MpSite("state-Im", nstore, input.prefix));
  StorageInitialize(sitesIm);

  DArray<2> a_subspace;
  DArray<2> b_subspace;
  DArray<2> s_subspace;
  DArray<2> d_subspace;

  vector<double> eigvs(2*nroots-1, 0.0);
  vector<double> eigvc(2*nroots-1, 0.0);
  vector<double> rnorm;
  vector<double> bnorm;

  double noise = input.noise;

  eigvs[0] = input.energy;
  eigvc[0] = input.energy;
  bool forward = 1;
  int nconv_roots = 1;
  int mroots = 1;

  //
  // linear-response
  //

  // TDA

  for(int iter = 0; iter < 1000; ++iter) {
    if(iter == 0) {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tMAKING GUESS WAVEFUNCTIONS FROM KRYLOV SUBSPACE "                                                   << endl;
      cout << "\t====================================================================================================" << endl;
    }
    else {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tDAVIDSON SWEEP ITERATION [ " << setw(4) << iter << " ] FOR TDA"   << endl;
      cout << "\t====================================================================================================" << endl;
    }

    mroots = lrt_sweep(1, sites0, sitesRe, sitesIm, a_subspace, b_subspace, s_subspace, d_subspace, eigvs, rnorm, bnorm, nroots, mroots, nconv_roots, TDA, noise);
//  lrt_sweep(forward, sites0, sitesRe, sitesIm, hsub, ssub, eigvs, rnorm, bnorm, nroots, mroots, nconv_roots);
//  forward = !forward;

    if(iter >  0) {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tDAVIDSON SWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
      cout.precision(16);
      cout << "\t\t\tEnergy[" << setw(2) << 0 << "] = " << setw(24) << fixed << eigvs[0] << endl;
      for(int i = 1; i < nroots; ++i) {
        cout.precision(16);
        cout << "\t\t\tEnergy[" << setw(2) << i << "] = " << setw(24) << fixed << eigvs[i] << " ( R-norm = ";
        cout.precision(2);
        cout << setw(8) << scientific << rnorm[i] << " / B-norm = " << setw(8) << scientific << bnorm[i] << " ) " << endl;
      }
      cout << "\t====================================================================================================" << endl;
      cout << endl;

      nconv_roots = 1;
      for(; nconv_roots < nroots; ++nconv_roots)
        if(fabs(rnorm[nconv_roots]) >= input.tolerance) break;
    }

    if(nconv_roots == nroots) break;

//  mroots += nroots - nconv_roots;
    if(mroots > mxritz) mroots = nroots;

// FOR TEST (don't perform correction)
//  mroots = nroots;
//  nconv_roots = nroots;
  }

  // RPA

  mroots = nroots;
  nconv_roots = 1;

  for(int iter = 0; iter < 1000; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tDAVIDSON SWEEP ITERATION [ " << setw(4) << iter << " ] FOR RPA"   << endl;
    cout << "\t====================================================================================================" << endl;

    mroots = lrt_sweep(1, sites0, sitesRe, sitesIm, a_subspace, b_subspace, s_subspace, d_subspace, eigvc, rnorm, bnorm, nroots, mroots, nconv_roots, RPA, noise);

    if(iter >  0) {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tDAVIDSON SWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
      cout.precision(16);
      cout << "\t\t\tEnergy[" << setw(2) << 0 << "] = " << setw(24) << fixed << eigvc[0] << endl;
      for(int i = 1; i < nroots; ++i) {
        cout.precision(16);
        cout << "\t\t\tEnergy[" << setw(2) << i << "] = " << setw(24) << fixed << eigvc[i] << " ( R-norm = ";
        cout.precision(2);
        cout << setw(8) << scientific << rnorm[i] << " / B-norm = " << setw(8) << scientific << bnorm[i] << " ) " << endl;
      }
      cout << "\t====================================================================================================" << endl;
      cout << endl;

      nconv_roots = 1;
      for(; nconv_roots < nroots; ++nconv_roots)
        if(fabs(rnorm[nconv_roots]) >= input.tolerance) break;
    }

    if(nconv_roots == nroots) break;

//  mroots += nroots - nconv_roots;
    if(mroots > mxritz) mroots = nroots;

// FOR TEST (don't perform correction)
//  mroots = nroots;
//  nconv_roots = nroots;
  }

  double ecorr1 = 0.0;
  double ecorr2 = 0.0;
  for(int iroot = 1; iroot < nroots; ++iroot) {
    ecorr1 -= eigvc[iroot] * bnorm[iroot];
    ecorr2 += 0.5 * (eigvc[iroot] - eigvs[iroot]);
  }

  cout.precision(16);
  cout << "\t====================================================================================================" << endl;
  for(int i = 1; i < nroots; ++i)
    cout << "\tTDA Excitation Energy[ " << setw(2) << i << " ] = " << setw(24) << fixed << eigvs[i] << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  for(int i = 1; i < nroots; ++i)
    cout << "\tRPA Excitation Energy[ " << setw(2) << i << " ] = " << setw(24) << fixed << eigvc[i] << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\tCorrelation Energy (Scheme 1) = " << setw(24) << fixed << ecorr1 << " ( used ) " << endl;
  cout << "\tCorrelation Energy (Scheme 2) = " << setw(24) << fixed << ecorr2 << endl;
  double ecorr = ecorr1;
  cout << "\t====================================================================================================" << endl;
  cout << "\t\tPRINTING RPA ENERGY" << endl;
  cout.precision(16);
  cout << "\t\t\tEnergy[  0 ] = " << setw(24) << fixed << eigvc[0] + ecorr << endl;
  for(int i = 1; i < nroots; ++i) {
    cout.precision(16);
    cout << "\t\t\tEnergy[ " << setw(2) << i << " ] = " << setw(24) << fixed << eigvc[i] + eigvc[0] + ecorr << endl;
  }
  cout << "\t====================================================================================================" << endl;

  return eigvc;
}



void prototype::compute_derived_storage
(bool forward, const DArray<4>& mpo, const DArray<3>& opr0, const DArray<3>& bra0, const DArray<3>& ket1, DArray<3>& fopr)
{
//cout << "debug: compute_derived_storage: opr0.shape() = " << opr0.shape() << endl;
//cout << "debug: compute_derived_storage: bra0.shape() = " << bra0.shape() << endl;
//cout << "debug: compute_derived_storage: ket1.shape() = " << ket1.shape() << endl;
  if(ket1.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, opr0, shape(0), bra0, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(0, 2), mpo, shape(0, 1), 1.0, scr2);
    fopr.free();
    Dcontract(1.0, scr2, shape(0, 2), ket1, shape(0, 1), 1.0, fopr);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, opr0, shape(0), bra0, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 0), mpo, shape(1, 3), 1.0, scr2);
    fopr.free();
    Dcontract(1.0, scr2, shape(3, 0), ket1, shape(1, 2), 1.0, fopr);
  }
}

void prototype::transfer_derived_storage
(bool forward, const DArray<4>& mpo, const DArray<3>& fopr, const DArray<3>& bra0, const DArray<3>& ket0, DArray<5>& gopr)
{
//cout << "debug: transfer_derived_storage: fopr.shape() = " << fopr.shape() << endl;
//cout << "debug: transfer_derived_storage: bra0.shape() = " << bra0.shape() << endl;
//cout << "debug: transfer_derived_storage: ket0.shape() = " << ket0.shape() << endl;
  if(fopr.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, fopr, shape(0), bra0, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(0, 2), mpo, shape(0, 1), 1.0, scr2);
    DArray<5> gtmp;
    Dcontract(1.0, scr2, shape(2), ket0, shape(1), 1.0, gtmp);
    Dpermute(gtmp, shape(0, 3, 1, 2, 4), gopr);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, fopr, shape(0), bra0, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 0), mpo, shape(1, 3), 1.0, scr2);
    DArray<5> gtmp;
    Dcontract(1.0, scr2, shape(3), ket0, shape(1), 1.0, gtmp);
    Dpermute(gtmp, shape(0, 4, 1, 2, 3), gopr);
  }
}

void prototype::transfer_derived_storage
(bool forward, const DArray<4>& mpo, const DArray<5>& fopr, const DArray<3>& bra0, const DArray<3>& ket0, DArray<5>& gopr)
{
//cout << "debug: transfer_derived_storage: fopr.shape() = " << fopr.shape() << endl;
//cout << "debug: transfer_derived_storage: bra0.shape() = " << bra0.shape() << endl;
//cout << "debug: transfer_derived_storage: ket0.shape() = " << ket0.shape() << endl;
  if(fopr.size() == 0) return;

  if(forward) {
    DArray<6> scr1;
    Dcontract(1.0, fopr, shape(2), bra0, shape(0), 1.0, scr1);
    DArray<6> scr2;
    Dcontract(1.0, scr1, shape(2, 4), mpo, shape(0, 1), 1.0, scr2);
    gopr.free();
    Dcontract(1.0, scr2, shape(2, 4), ket0, shape(0, 1), 1.0, gopr);
  }
  else {
    DArray<6> scr1;
    Dcontract(1.0, fopr, shape(2), bra0, shape(2), 1.0, scr1);
    DArray<6> scr2;
    Dcontract(1.0, scr1, shape(5, 2), mpo, shape(1, 3), 1.0, scr2);
    gopr.free();
    Dcontract(1.0, scr2, shape(5, 2), ket0, shape(1, 2), 1.0, gopr);
  }
}

void prototype::compute_a_diagonal_elements
(bool forward, double e0, const DArray<4>& mpo, const DArray<3>& opr0, const DArray<3>& copr,
               const DArray<3>& bra1, const DArray<3>& ket1, DArray<4>& aii)
{
//cout << "debug: compute_a_diagonal: opr0.shape() = " << opr0.shape() << endl;
//cout << "debug: compute_a_diagonal: copr.shape() = " << copr.shape() << endl;
//cout << "debug: compute_a_diagonal: bra1.shape() = " << bra1.shape() << endl;
//cout << "debug: compute_a_diagonal: ket1.shape() = " << ket1.shape() << endl;
  if(bra1.size() == 0 || ket1.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, opr0, shape(0), bra1, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(0, 2), mpo, shape(0, 1), 1.0, scr2);
    DArray<3> scr3;
    Dcontract(1.0, scr2, shape(0, 2), ket1, shape(0, 1), 1.0, scr3);
    DArray<4> aiitmp;
    Dcontract(1.0, scr3, shape(1), copr, shape(1), 1.0, aiitmp);
    Dpermute(aiitmp, shape(0, 2, 1, 3), aii);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, opr0, shape(0), bra1, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 0), mpo, shape(1, 3), 1.0, scr2);
    DArray<3> scr3;
    Dcontract(1.0, scr2, shape(3, 0), ket1, shape(1, 2), 1.0, scr3);
    DArray<4> aiitmp;
    Dcontract(1.0, scr3, shape(1), copr, shape(1), 1.0, aiitmp);
    Dpermute(aiitmp, shape(0, 2, 1, 3), aii);
  }
  for(int i0 = 0; i0 < aii.extent(0); ++i0) {
    for(int i1 = 0; i1 < aii.extent(1); ++i1) {
      aii(i0, i1, i0, i1) -= e0;
    }
  }
}

void prototype::compute_a_nearest_elements
(bool forward, const DArray<4>& mpo, const DArray<3>& fopr, const DArray<3>& copr,
               const DArray<3>& bra1, const DArray<3>& ket0, DArray<4>& aij)
{
//cout << "debug: compute_a_nearest: fopr.shape() = " << fopr.shape() << endl;
//cout << "debug: compute_a_nearest: copr.shape() = " << copr.shape() << endl;
//cout << "debug: compute_a_nearest: bra1.shape() = " << bra1.shape() << endl;
//cout << "debug: compute_a_nearest: ket0.shape() = " << ket0.shape() << endl;
  if(fopr.size() == 0 || bra1.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(2), ket0, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 1), mpo, shape(2, 3), 1.0, scr2);
    DArray<4> scr3;
    Dcontract(1.0, fopr, shape(0), bra1, shape(0), 1.0, scr3);
    DArray<4> aijtmp;
    Dcontract(1.0, scr3, shape(0, 2), scr2, shape(2, 3), 1.0, aijtmp);
    Dpermute(aijtmp, shape(0, 3, 1, 2), aij);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(2), ket0, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(1, 2), mpo, shape(0, 2), 1.0, scr2);
    DArray<4> scr3;
    Dcontract(1.0, fopr, shape(0), bra1, shape(2), 1.0, scr3);
    DArray<4> aijtmp;
    Dcontract(1.0, scr3, shape(3, 0), scr2, shape(2, 3), 1.0, aijtmp);
    Dpermute(aijtmp, shape(0, 3, 1, 2), aij);
  }
}

void prototype::compute_a_matrix_elements
(bool forward, const DArray<4>& mpo, const DArray<5>& fopr, const DArray<3>& copr,
               const DArray<3>& bra1, const DArray<3>& ket0, DArray<4>& aij)
{
//cout << "debug: compute_a_matrix: fopr.shape() = " << fopr.shape() << endl;
//cout << "debug: compute_a_matrix: copr.shape() = " << copr.shape() << endl;
//cout << "debug: compute_a_matrix: bra1.shape() = " << bra1.shape() << endl;
//cout << "debug: compute_a_matrix: ket0.shape() = " << ket0.shape() << endl;
  if(fopr.size() == 0 || bra1.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(2), ket0, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 1), mpo, shape(2, 3), 1.0, scr2);
    DArray<6> scr3;
    Dcontract(1.0, fopr, shape(2), bra1, shape(0), 1.0, scr3);
    aij.free();
    Dcontract(1.0, scr3, shape(4, 2, 3), scr2, shape(3, 2, 1), 1.0, aij);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(2), ket0, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(1, 2), mpo, shape(0, 2), 1.0, scr2);
    DArray<6> scr3;
    Dcontract(1.0, fopr, shape(2), bra1, shape(2), 1.0, scr3);
    aij.free();
    Dcontract(1.0, scr3, shape(5, 2, 3), scr2, shape(2, 3, 1), 1.0, aij);
  }
}

void prototype::compute_b_nearest_elements
(bool forward, const DArray<4>& mpo, const DArray<3>& fopr, const DArray<3>& copr,
               const DArray<3>& bra0, const DArray<3>& ket1, DArray<4>& bij)
{
  if(fopr.size() == 0 || ket1.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(0), bra0, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 0), mpo, shape(1, 3), 1.0, scr2);
    DArray<3> scr3;
    Dcontract(1.0, fopr, shape(0, 1), scr2, shape(1, 2), 1.0, scr3);
    DArray<4> bijtmp;
    Dcontract(1.0, scr3, shape(2), ket1, shape(1), 1.0, bijtmp);
    Dpermute(bijtmp, shape(0, 2, 3, 1), bij);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(0), bra0, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(0, 2), mpo, shape(0, 1), 1.0, scr2);
    DArray<3> scr3;
    Dcontract(1.0, fopr, shape(0, 1), scr2, shape(1, 3), 1.0, scr3);
    DArray<4> bijtmp;
    Dcontract(1.0, scr3, shape(2), ket1, shape(1), 1.0, bijtmp);
    Dpermute(bijtmp, shape(0, 3, 2, 1), bij);
  }
}

void prototype::compute_b_matrix_elements
(bool forward, const DArray<4>& mpo, const DArray<5>& fopr, const DArray<3>& copr,
               const DArray<3>& bra0, const DArray<3>& ket1, DArray<4>& bij)
{
  if(fopr.size() == 0 || ket1.size() == 0) return;

  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(0), bra0, shape(2), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(3, 0), mpo, shape(1, 3), 1.0, scr2);
    DArray<6> scr3;
    Dcontract(1.0, fopr, shape(4), ket1, shape(0), 1.0, scr3);
    bij.free();
    Dcontract(1.0, scr3, shape(2, 3, 4), scr2, shape(1, 2, 3), 1.0, bij);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, copr, shape(0), bra0, shape(0), 1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(0, 2), mpo, shape(0, 1), 1.0, scr2);
    DArray<6> scr3;
    Dcontract(1.0, fopr, shape(4), ket1, shape(2), 1.0, scr3);
    bij.free();
    Dcontract(1.0, scr3, shape(2, 3, 5), scr2, shape(1, 3, 2), 1.0, bij);
  }
}

void prototype::compute_s_matrix_elements
(bool forward, const DArray<2>& sopr, const DArray<3>& bra1, const DArray<3>& ket1, DArray<4>& sii)
{
  if(bra1.size() == 0 || ket1.size() == 0) return;

  if(forward) {
    DArray<2> scr1;
    Dcontract(1.0, bra1, shape(0, 1), ket1, shape(0, 1), 1.0, scr1);
    DArray<4> siitmp;
    Dger(1.0, scr1, sopr, siitmp);
    Dpermute(siitmp, shape(0, 2, 1, 3), sii);
  }
  else {
    DArray<2> scr1;
    Dcontract(1.0, bra1, shape(1, 2), ket1, shape(1, 2), 1.0, scr1);
    DArray<4> siitmp;
    Dger(1.0, scr1, sopr, siitmp);
    Dpermute(siitmp, shape(0, 2, 1, 3), sii);
  }
}

void prototype::construct_effective_hamiltonian(MpStorages& sites, double eigenvalue)
{
  int L = sites.size();
  vector< DArray<4> > a_data(L*L, DArray<4>());
  vector< DArray<4> > b_data(L*L, DArray<4>());
  vector< DArray<4> > s_data(L,   DArray<4>());

  sites[0].load(0);
  for(int i = 0; i < L-1; ++i) {
    sites[i+1].load(i+1);
    Canonicalize(1, sites[i].wfnc[0], sites[i].lmps[0], sites[i].ltnj[0]);
    sites[i+1].wfnc[0].free();
    ComputeGuess(1, sites[i].lmps[0], sites[i].wfnc[0], sites[i+1].rmps[0], sites[i+1].wfnc[0]);
    sites[i+1].lopr[0].free();
    Renormalize(1, sites[i].mpo,  sites[i].lopr[0], sites[i].lmps[0], sites[i].lmps[0], sites[i+1].lopr[0]);
    sites[i].save(i);
  }
  Canonicalize(1, sites[L-1].wfnc[0], sites[L-1].lmps[0], sites[L-1].ltnj[0]);

  vector< DArray<3> > copr(L, DArray<3>());
  copr[L-1].resize(1, 1, 1); copr[L-1] = 1.0;

  vector< DArray<2> > sopr(L, DArray<2>());
  sopr[L-1].resize(1, 1); sopr[L-1] = 1.0;

  Dcopy(sites[L-1].wfnc[0], sites[L-1].lmps[0]);
  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    Renormalize(0, sites[i].mpo, copr[i], sites[i].lmps[0], sites[i].lmps[0], copr[i-1]);
    Renormalize(0,               sopr[i], sites[i].lmps[0], sites[i].lmps[0], sopr[i-1]);
    sites[i].save(i);
  }

  vector< DArray<3> > fopr(L, DArray<3>());
  vector< DArray<5> > gopr(L, DArray<5>());

  vector<int> nrows(L, 0);
  vector<int> ncols(L, 0);

  for(int i = 0; i < L; ++i) {
    cout << "\tcomputing matrix elements @ site[" << i << "]" << endl;
    if(i > 0) sites[i].load(i);

    int ii = i*L+i;
    compute_a_diagonal_elements(1, eigenvalue, sites[i].mpo, sites[i].lopr[0], copr[i], sites[i].ltnj[0], sites[i].ltnj[0], a_data[ii]);
    compute_s_matrix_elements(1, sopr[i], sites[i].ltnj[0], sites[i].ltnj[0], s_data[i]);

    nrows[i] = a_data[ii].extent(0)*a_data[ii].extent(1);
    ncols[i] = a_data[ii].extent(2)*a_data[ii].extent(3);

    fopr[i].free();
    compute_derived_storage(1, sites[i].mpo, sites[i].lopr[0], sites[i].lmps[0], sites[i].ltnj[0], fopr[i]);

    for(int j = 0; j < i-1; ++j) {
      cout << "\tcomputing matrix elements @ site[" << i << "] - site[" << j << "]" << endl;
      int ij = i*L+j;
      int ji = j*L+i;

      compute_a_matrix_elements(1, sites[i].mpo, gopr[j], copr[i], sites[i].ltnj[0], sites[i].lmps[0], a_data[ij]);
      if(a_data[ij].size() > 0) Dpermute(a_data[ij], shape(2, 3, 0, 1), a_data[ji]);

      compute_b_matrix_elements(1, sites[i].mpo, gopr[j], copr[i], sites[i].lmps[0], sites[i].ltnj[0], b_data[ij]);
      if(b_data[ij].size() > 0) Dpermute(b_data[ij], shape(2, 3, 0, 1), b_data[ji]);

      DArray<5> gscr;
      transfer_derived_storage(1, sites[i].mpo, gopr[j], sites[i].lmps[0], sites[i].lmps[0], gscr);
      Dcopy(gscr, gopr[j]);
    }

    if(i > 0) {
      int j = i-1;

      cout << "\tcomputing matrix elements @ site[" << i << "] - site[" << j << "]" << endl;
      int ij = i*L+j;
      int ji = j*L+i;

      compute_a_nearest_elements(1, sites[i].mpo, fopr[j], copr[i], sites[i].ltnj[0], sites[i].lmps[0], a_data[ij]);
      if(a_data[ij].size() > 0) Dpermute(a_data[ij], shape(2, 3, 0, 1), a_data[ji]);

      compute_b_nearest_elements(1, sites[i].mpo, fopr[j], copr[i], sites[i].lmps[0], sites[i].ltnj[0], b_data[ij]);
      if(b_data[ij].size() > 0) Dpermute(b_data[ij], shape(2, 3, 0, 1), b_data[ji]);

      DArray<5> gscr;
      transfer_derived_storage(1, sites[i].mpo, fopr[j], sites[i].lmps[0], sites[i].lmps[0], gscr);
      Dcopy(gscr, gopr[j]);
    }

    sites[i].save(i);
  }

  int nrows_full = accumulate(nrows.begin(), nrows.end(), 0);
  int ncols_full = accumulate(ncols.begin(), ncols.end(), 0);

  DArray<2> full_a_matrix(nrows_full, ncols_full); full_a_matrix = 0.0;
  DArray<2> full_b_matrix(nrows_full, ncols_full); full_b_matrix = 0.0;
  DArray<2> full_s_matrix(nrows_full, ncols_full); full_s_matrix = 0.0;
  DArray<2> full_d_matrix(nrows_full, ncols_full); full_d_matrix = 0.0;

  int ibeg = 0;
  for(int i = 0; i < L; ++i) {
    int iend = ibeg+nrows[i]-1;
    int jbeg = 0;
    for(int j = 0; j < L; ++j) {
      int jend = jbeg+ncols[j]-1;
      int ij = i*L+j;
      if(a_data[ij].size() > 0)
        Dcopy_direct(a_data[ij], full_a_matrix(Range(ibeg, iend), Range(jbeg, jend)));
      if(b_data[ij].size() > 0)
        Dcopy_direct(b_data[ij], full_b_matrix(Range(ibeg, iend), Range(jbeg, jend)));
      jbeg += ncols[j];
    }
    if(s_data[i].size() > 0)
      Dcopy_direct(s_data[i], full_s_matrix(Range(ibeg, iend), Range(ibeg, iend)));

    ibeg += nrows[i];
  }

  Dscal(1.0/L, full_a_matrix);
  Dscal(1.0/L, full_b_matrix);
  Dscal(1.0/L, full_s_matrix);

  cout << "\t********** printing A-matrix **********" << endl;
  printing_symmetric_matrix(full_a_matrix);
  cout << "\t********** printing B-matrix **********" << endl;
  printing_symmetric_matrix(full_b_matrix);
  cout << "\t********** printing S-matrix **********" << endl;
  printing_symmetric_matrix(full_s_matrix);

  int nrows_ab = full_a_matrix.rows() + full_b_matrix.rows();
  int ncols_ab = full_a_matrix.cols() + full_b_matrix.cols();
  DArray<2> st_matrix(nrows_ab, ncols_ab); st_matrix = 0.0; // stability matrix
  DArray<2> mt_matrix(nrows_ab, ncols_ab); mt_matrix = 0.0; // metric
  for(int i = 0; i < nrows_full; ++i) {
    for(int j = 0; j < ncols_full; ++j) {
      st_matrix(i,            j           ) = full_a_matrix(i, j);
      st_matrix(i+nrows_full, j           ) = full_b_matrix(j, i);
      st_matrix(i,            j+ncols_full) = full_b_matrix(i, j);
      st_matrix(i+nrows_full, j+ncols_full) = full_a_matrix(j, i);

      mt_matrix(i,            j           ) = full_s_matrix(i, j);
      mt_matrix(i+nrows_full, j+ncols_full) =-full_s_matrix(j, i);
    }
  }

  {

  DArray<1> freq;
  DArray<2> alpha;
  ComputeEigenvalues(full_a_matrix, full_s_matrix, freq, alpha);
  cout << "\tprinting TDA frequencies" << endl;
  for(int i = 0; i < freq.size(); ++i)
    cout << "\t" << setw(12) << setprecision(8) << fixed << freq(i) << endl;

  }

  {

  DArray<1> freq;
  DArray<2> alphaRe;
  DArray<2> alphaIm;
  ComputeEigenvalues(full_a_matrix, full_b_matrix, full_s_matrix, full_d_matrix, freq, alphaRe, alphaIm);
  cout << "\tprinting RPA frequencies" << endl;
  for(int i = 0; i < freq.size(); ++i)
    cout << "\t" << setw(12) << setprecision(8) << fixed << freq(i) << endl;

  }
}

void prototype::printing_symmetric_matrix(const DArray<2>& a, int mcols)
{
  int nrows = a.rows();
  int ncols = a.cols();
  int ndivs = ncols / mcols;
  int nmods = ncols % mcols;
  for(int idivs = 0; idivs < ndivs; ++idivs) {
    int jbeg = idivs*mcols;
    int jend = jbeg +mcols;
    cout << "------------" << "\t";
    for(int j = jbeg; j < jend; ++j) {
      cout << setw(12) << j;
    }
    cout << endl << endl;
    for(int i = jbeg; i < nrows; ++i) {
      cout << setw(12) << i << "\t";
      int jtrm = min(i+1, jend);
      for(int j = jbeg; j < jtrm; ++j) {
        cout << setw(12) << setprecision(8) << fixed << a(i, j);
      }
      cout << endl;
    }
    cout << endl;
  }
  if(nmods == 0) return;
  {
    int jbeg = ndivs*mcols;
    int jend = jbeg +nmods;
    cout << "------------" << "\t";
    for(int j = jbeg; j < jend; ++j) {
      cout << setw(12) << j;
    }
    cout << endl << endl;
    for(int i = jbeg; i < nrows; ++i) {
      cout << setw(12) << i << "\t";
      int jtrm = min(i+1, jend);
      for(int j = jbeg; j < jtrm; ++j) {
        cout << setw(12) << setprecision(8) << fixed << a(i, j);
      }
      cout << endl;
    }
    cout << endl;
  }
}








