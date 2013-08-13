#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; };

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
  for(int i = 0; i < h_subspace.shape(0); ++i) {
    cout << "\t";
    for(int j = 0; j < h_subspace.shape(1); ++j) {
      cout << setw(12) << fixed << h_subspace(i, j);
    }
    cout << endl;
  }
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\t\t\tPRINTING SUBSPACE OVERLAP" << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout.precision(8);
  for(int i = 0; i < s_subspace.shape(0); ++i) {
    cout << "\t";
    for(int j = 0; j < s_subspace.shape(1); ++j) {
      cout << setw(12) << fixed << s_subspace(i, j);
    }
    cout << endl;
  }
}

void printing_matrix_elements(const DArray<2>& a_subspace, const DArray<2>& b_subspace,
                              const DArray<2>& s_subspace, const DArray<2>& d_subspace)
{
  int nrows = a_subspace.shape(0);
  int ncols = a_subspace.shape(1);
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
    QSDArray<4> LHS; // LHS = C(I)R(0) + C(0)R(I)
    QSDgemm(NoTrans, NoTrans, 1.0, lsite1.wfnc[iroot], rsite0.rmps[0], 1.0, LHS);
    QSDgemm(NoTrans, NoTrans, 1.0, lsite0.wfnc[0], rsite1.rmps[iroot], 1.0, LHS);

    QSDArray<4> RHS; // RHS = L(I)C(0) + L(0)C(I)
    QSDgemm(NoTrans, NoTrans, 1.0, lsite1.lmps[iroot], rsite0.wfnc[0], 1.0, RHS);
    QSDgemm(NoTrans, NoTrans, 1.0, lsite0.lmps[0], rsite1.wfnc[iroot], 1.0, RHS);

    QSDaxpy(-1.0, RHS, LHS); // LHS = LHS - RHS
    cout << "\tChecking gauge-condition[" << setw(2) << iroot << "]: subnorm = " << setw(20) << fixed << QSDdotc(LHS, LHS) << endl;
  }

  cout << "\tChecking site-exchange: " << endl;
  {
    QSDArray<4> psi;
    QSDgemm(NoTrans, NoTrans, 1.0, lsite0.wfnc[0], rsite0.rmps[0], 1.0, psi);
    {
       SDArray<1> s;
      QSDArray<3> u;
      QSDArray<3> v;
      QSDgesvd(LeftArrow, psi, s, u, v, 0);
      cout << "\t\tSchmidt values (original): " << s << endl;
    }
    {
      QSDArray<4> psi_pmt;
      QSDpermute(psi, shape(0, 2, 1, 3), psi_pmt);
       SDArray<1> s;
      QSDArray<3> u;
      QSDArray<3> v;
      QSDgesvd(LeftArrow, psi_pmt, s, u, v, 0);
      cout << "\t\tSchmidt values (exchange): " << s << endl;
    }
  }
}

void prototype::RotateStorage
(const DArray<2>& alpha, vector< QSDArray<3> >& store, int nroots)
{
  int mroots = 1 + alpha.shape(0);
  assert(store.size() >= mroots);
  vector< QSDArray<3> > tmpstr(nroots);
  for(int iroot = 1; iroot < nroots; ++iroot)
    for(int jroot = 1; jroot < mroots; ++jroot)
      QSDaxpy(alpha(jroot-1, iroot-1), store[jroot], tmpstr[iroot]);

  for(int iroot = 1; iroot < nroots; ++iroot)
    QSDcopy(tmpstr[iroot], store[iroot]);
}

void prototype::RotateStorage
(const DArray<2>& alphaRe, const DArray<2>& alphaIm, vector< QSDArray<3> >& storeRe, vector< QSDArray<3> >& storeIm, int nroots)
{
  int mroots = 1 + alphaRe.shape(0);
  assert(storeRe.size() >= mroots);
  assert(storeIm.size() >= mroots);
  vector< QSDArray<3> > tmpstrRe(nroots);
  vector< QSDArray<3> > tmpstrIm(nroots);
  for(int iroot = 1; iroot < nroots; ++iroot) {
    for(int jroot = 1; jroot < mroots; ++jroot) {
      QSDaxpy(alphaRe(jroot-1, iroot-1), storeRe[jroot], tmpstrRe[iroot]);
      QSDaxpy(alphaIm(jroot-1, iroot-1), storeIm[jroot], tmpstrRe[iroot]);
      QSDaxpy(alphaRe(jroot-1, iroot-1), storeIm[jroot], tmpstrIm[iroot]);
      QSDaxpy(alphaIm(jroot-1, iroot-1), storeRe[jroot], tmpstrIm[iroot]);
    }
  }

  for(int iroot = 1; iroot < nroots; ++iroot) {
    QSDcopy(tmpstrRe[iroot], storeRe[iroot]);
    QSDcopy(tmpstrIm[iroot], storeIm[iroot]);
  }
}

void prototype::SolveCorrectionEquation
(bool forward, const MpSite& site0, MpSite& siteRe, MpSite& siteIm, const QSDArray<3>& diag,
 const vector<double>& eigv, vector<double>& rnorm, int nroots, int mroots, int kroots, bool boundary, const CALC_TYPE& calc_type)
{
  if(mroots == 1) {
    QSDArray<3> tmpwfc(site0.wfnc[0]);
    for(int iroot = 1; iroot < nroots; ++iroot) {
      siteRe.wfnc[iroot].clear();
      ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], tmpwfc, siteRe.wfnc[iroot]);

      Orthogonalize(site0.wfnc[0], siteRe.wfnc[iroot]);
      for(int jroot = 1; jroot < iroot; ++jroot) {
        Normalize(siteRe.wfnc[iroot]);
        Orthogonalize(site0.wfnc[0], siteRe.wfnc[iroot]);
      }
      Normalize(siteRe.wfnc[iroot]);
      QSDcopy(siteRe.wfnc[iroot], tmpwfc);

      QSDcopy(site0.wfnc[0], siteIm.wfnc[iroot]);
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

  vector< QSDArray<3> > sgvxRe(nroots);
  vector< QSDArray<3> > sgvxIm(nroots);

  // diagonal block
  for(int iroot = 1; iroot < nroots; ++iroot) {
    ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], siteRe.wfnc[iroot], sgvxRe[iroot]);
    QSDaxpy(-eigv[0], siteRe.wfnc[iroot], sgvxRe[iroot]);
    ComputeSigmaVector(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxRe[iroot]);
    ComputeSigmaVector(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgvxRe[iroot]);
  }
  for(int iroot = 1; iroot < nroots; ++iroot) {
    ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], siteIm.wfnc[iroot], sgvxIm[iroot]);
    QSDaxpy(-eigv[0], siteIm.wfnc[iroot], sgvxIm[iroot]);
    ComputeSigmaVector(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxIm[iroot]);
    ComputeSigmaVector(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgvxIm[iroot]);
  }
  // off-diagonal block
  if(calc_type == RPA) {

  for(int iroot = 1; iroot < nroots; ++iroot) {
    ComputeSigmaVectorConj(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxRe[iroot]);
    ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgvxRe[iroot]);
  }
  for(int iroot = 1; iroot < nroots; ++iroot) {
    ComputeSigmaVectorConj(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgvxIm[iroot]);
    ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgvxIm[iroot]);
  }

  }

  // residual norm & correction
  for(int iroot = 1; iroot < nroots; ++iroot) {
    QSDArray<3> errvRe(sgvxRe[iroot]);
    QSDaxpy(-eigv[iroot], siteRe.wfnc[iroot], errvRe);

    QSDArray<3> errvIm(sgvxIm[iroot]);
    QSDaxpy(+eigv[iroot], siteIm.wfnc[iroot], errvIm);

    rnorm[iroot] += (QSDdotc(errvRe, errvRe) + QSDdotc(errvIm, errvIm));
//  rnorm[iroot] += QSDdotc(errvRe, errvRe);

    if(iroot < kroots) continue;

    davidson::precondition(+eigv[iroot], diag, errvRe);
    davidson::precondition(-eigv[iroot], diag, errvIm);

//  if(!boundary) {
    if(0) {
      if(forward) {
        Orthogonalize(0, site0.rmps[0], errvRe);
        Orthogonalize(0, site0.rmps[0], errvIm);
      }
      else {
        Orthogonalize(1, site0.lmps[0], errvRe);
        Orthogonalize(1, site0.lmps[0], errvIm);
      }
      QSDaxpy(1.0, siteRe.wfnc[mroots+iroot-kroots], errvRe);
      QSDaxpy(1.0, siteIm.wfnc[mroots+iroot-kroots], errvIm);
    }

    Orthogonalize(site0.wfnc[0], errvRe);
    Orthogonalize(site0.wfnc[0], errvIm);

    double normRe = QSDdotc(errvRe, errvRe);
    double normIm = QSDdotc(errvIm, errvIm);

//  if(normRe - normIm < 1.0e-8) {
//    cout << "\t\tIgnore correction for imaginary wavefunction, since R-norm is larger than real part" << endl;
//    errvIm = 0.0;
//  }
//  if(normRe > 1.0) {
      QSDscal(1.0/sqrt(normRe), errvRe);
      QSDscal(1.0/sqrt(normRe), errvIm);
//  }

//  if(normRe - normIm < 0.0) {
//    QSDscal(-1.0, errvIm);
//    QSDcopy(errvIm, siteRe.wfnc[mroots+iroot-kroots]);
//    QSDcopy(errvRe, siteIm.wfnc[mroots+iroot-kroots]);
//  }
//  else {
      QSDcopy(errvRe, siteRe.wfnc[mroots+iroot-kroots]);
      QSDcopy(errvIm, siteIm.wfnc[mroots+iroot-kroots]);
//  }
  }
}

void prototype::ComputeMatrixElements
(bool forward, const MpSite& site0, const MpSite& siteRe, const MpSite& siteIm,
 const vector<double>& eigv, DArray<2>& a_subspace, DArray<2>& b_subspace, DArray<2>& s_subspace, DArray<2>& d_subspace,
       vector<double>& bnorm, int nroots, int mroots, int kroots, bool boundary, bool last, const CALC_TYPE& calc_type)
{
  int lroots = mroots + nroots - kroots;

  // projection to orthogonal subspace
  vector< QSDArray<3> > qfncRe(lroots);
  vector< QSDArray<3> > qfncIm(lroots);
  for(int iroot = 1; iroot < lroots; ++iroot) {
    QSDcopy(siteRe.wfnc[iroot], qfncRe[iroot]);
    QSDcopy(siteIm.wfnc[iroot], qfncIm[iroot]);

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
    double sii = QSDdotc(qfncRe[iroot], qfncRe[iroot]);
    s_subspace(iroot-1, iroot-1) += sii;
    for(int jroot = 1; jroot < iroot; ++jroot) {
      double sij = QSDdotc(qfncRe[iroot], qfncRe[jroot]);
      s_subspace(iroot-1, jroot-1) += sij;
      s_subspace(jroot-1, iroot-1) += sij;
    }
  }
  // Yi Sij Yj
  for(int iroot = 1; iroot < lroots; ++iroot) {
    double sii = QSDdotc(qfncIm[iroot], qfncIm[iroot]);
    bnorm[iroot] += sii; // save Y-norm to compute RPA energy
    s_subspace(iroot-1, iroot-1) -= sii;
    for(int jroot = 1; jroot < iroot; ++jroot) {
      double sij = QSDdotc(qfncIm[iroot], qfncIm[jroot]);
      s_subspace(iroot-1, jroot-1) -= sij;
      s_subspace(jroot-1, iroot-1) -= sij;
    }
  }
  // Xi Sij Yj
  for(int iroot = 1; iroot < lroots; ++iroot) {
    for(int jroot = 1; jroot < iroot; ++jroot) {
      double dij = QSDdotc(qfncRe[iroot], qfncIm[jroot]);
      double dji = QSDdotc(qfncRe[jroot], qfncIm[iroot]);
      d_subspace(iroot-1, jroot-1) += (dij - dji);
      d_subspace(jroot-1, iroot-1) += (dji - dij);
    }
  }
  // compute reduced A-matrix & B-matrix / 1-st order terms
  // Aij += Xi Hij Xj + Yi Hij Yj, Bij = Xi Hij Yj + Yi Hij Xj
  {
    vector< QSDArray<3> > sgvxRe(lroots);
    vector< QSDArray<3> > sgv0Re(lroots);

    vector< QSDArray<3> > sgvxIm(lroots);
    vector< QSDArray<3> > sgv0Im(lroots);

    for(int iroot = 1; iroot < lroots; ++iroot) {
        ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], qfncRe[iroot], sgvxRe[iroot]);
        QSDaxpy(-eigv[0], qfncRe[iroot], sgvxRe[iroot]);
      if(forward)
        ComputeSigmaVector(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Re[iroot]);
      else
        ComputeSigmaVector(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgv0Re[iroot]);
    }
    for(int iroot = 1; iroot < lroots; ++iroot) {
        ComputeSigmaVector(site0.mpo, site0.lopr[0], site0.ropr[0], qfncIm[iroot], sgvxIm[iroot]);
        QSDaxpy(-eigv[0], qfncIm[iroot], sgvxIm[iroot]);
      if(forward)
        ComputeSigmaVector(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Im[iroot]);
      else
        ComputeSigmaVector(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgv0Im[iroot]);
    }

//cout.precision(8);
//cout << "debug: check a_subspace (real x real)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = QSDdotc(qfncRe[iroot], sgvxRe[iroot])
                 + QSDdotc(qfncRe[iroot], sgv0Re[iroot]) * 2.0;
      a_subspace(iroot-1, iroot-1) += vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = QSDdotc(qfncRe[iroot], sgvxRe[jroot])
                   + QSDdotc(qfncRe[iroot], sgv0Re[jroot])
                   + QSDdotc(qfncRe[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        a_subspace(iroot-1, jroot-1) += vij;
        a_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << vii << endl;
    }
//cout << "debug: check a_subspace (imag x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = QSDdotc(qfncIm[iroot], sgvxIm[iroot])
                 + QSDdotc(qfncIm[iroot], sgv0Im[iroot]) * 2.0;
      a_subspace(iroot-1, iroot-1) += vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = QSDdotc(qfncIm[iroot], sgvxIm[jroot])
                   + QSDdotc(qfncIm[iroot], sgv0Im[jroot])
                   + QSDdotc(qfncIm[jroot], sgv0Im[iroot]);
//cout << setw(12) << vij;
        a_subspace(iroot-1, jroot-1) += vij;
        a_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << vii << endl;
    }
//cout << "debug: check b_subspace (real x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = QSDdotc(qfncRe[iroot], sgvxIm[iroot])
                 + QSDdotc(qfncRe[iroot], sgv0Im[iroot])
                 + QSDdotc(qfncIm[iroot], sgv0Re[iroot]);
      b_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = QSDdotc(qfncRe[iroot], sgvxIm[jroot])
                   + QSDdotc(qfncRe[iroot], sgv0Im[jroot])
                   + QSDdotc(qfncRe[jroot], sgv0Im[iroot])
                   + QSDdotc(qfncIm[iroot], sgvxRe[jroot])
                   + QSDdotc(qfncIm[iroot], sgv0Re[jroot])
                   + QSDdotc(qfncIm[jroot], sgv0Re[iroot]);
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

  QSDArray<2> gauge0;
  if(forward)
    QSDgemm(NoTrans, ConjTrans, 1.0, site0.wfnc[0], site0.rmps[0], 1.0, gauge0);
  else
    QSDgemm(ConjTrans, NoTrans, 1.0, site0.lmps[0], site0.wfnc[0], 1.0, gauge0);

  // compute reduced A-matrix & B-matrix / 2-nd order terms
  // Aij += Xi Wij Yj + Yi Wij Xj, Bij = Xi Wij Xj + Yi Wij Yj
  {
    vector< QSDArray<3> > sgv0Re(lroots);
    vector< QSDArray<3> > sgv0Im(lroots);

    for(int iroot = 1; iroot < lroots; ++iroot) {
      qfncRe[iroot].clear();
      if(forward) {
        QSDgemm(NoTrans, NoTrans, 1.0, gauge0, siteRe.rmps[iroot], 1.0, qfncRe[iroot]);
        Orthogonalize(site0.wfnc[0], qfncRe[iroot]);
        ComputeSigmaVectorConj(site0.mpo, siteRe.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Re[iroot]);
//cout << "debug: forward: qfncRe[" << iroot << "] = " << qfncRe[iroot] << endl;
      }
      else {
        QSDgemm(NoTrans, NoTrans, 1.0, siteRe.lmps[iroot], gauge0, 1.0, qfncRe[iroot]);
        Orthogonalize(site0.wfnc[0], qfncRe[iroot]);
        ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteRe.ropr[iroot], site0.wfnc[0], sgv0Re[iroot]);
//cout << "debug: backward: qfncRe[" << iroot << "] = " << qfncRe[iroot] << endl;
      }
    }
    for(int iroot = 1; iroot < lroots; ++iroot) {
      qfncIm[iroot].clear();
      if(forward) {
        QSDgemm(NoTrans, NoTrans, 1.0, gauge0, siteIm.rmps[iroot], 1.0, qfncIm[iroot]);
        Orthogonalize(site0.wfnc[0], qfncIm[iroot]);
        ComputeSigmaVectorConj(site0.mpo, siteIm.lopr[iroot], site0.ropr[0], site0.wfnc[0], sgv0Im[iroot]);
//cout << "debug: forward: qfncIm[" << iroot << "] = " << qfncIm[iroot] << endl;
      }
      else {
        QSDgemm(NoTrans, NoTrans, 1.0, siteIm.lmps[iroot], gauge0, 1.0, qfncIm[iroot]);
        Orthogonalize(site0.wfnc[0], qfncIm[iroot]);
        ComputeSigmaVectorConj(site0.mpo, site0.lopr[0], siteIm.ropr[iroot], site0.wfnc[0], sgv0Im[iroot]);
//cout << "debug: backward: qfncIm[" << iroot << "] = " << qfncIm[iroot] << endl;
      }
    }

//cout << "debug: check a_subspace (real x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = QSDdotc(qfncRe[iroot], sgv0Im[iroot])
                 + QSDdotc(qfncIm[iroot], sgv0Re[iroot]);
      a_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = QSDdotc(qfncRe[iroot], sgv0Im[jroot])
                   + QSDdotc(qfncRe[jroot], sgv0Im[iroot])
                   + QSDdotc(qfncIm[iroot], sgv0Re[jroot])
                   + QSDdotc(qfncIm[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        a_subspace(iroot-1, jroot-1) += vij;
        a_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << 2.0 * vii << endl;
    }
//cout << "debug: check b_subspace (real x real)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = QSDdotc(qfncRe[iroot], sgv0Re[iroot]);
      b_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = QSDdotc(qfncRe[iroot], sgv0Re[jroot])
                   + QSDdotc(qfncRe[jroot], sgv0Re[iroot]);
//cout << setw(12) << vij;
        b_subspace(iroot-1, jroot-1) += vij;
        b_subspace(jroot-1, iroot-1) += vij;
      }
//cout << setw(12) << 2.0 * vii << endl;
    }
//cout << "debug: check b_subspace (imag x imag)" << endl;
    for(int iroot = 1; iroot < lroots; ++iroot) {
      double vii = QSDdotc(qfncIm[iroot], sgv0Im[iroot]);
      b_subspace(iroot-1, iroot-1) += 2.0 * vii;
      for(int jroot = 1; jroot < iroot; ++jroot) {
        double vij = QSDdotc(qfncIm[iroot], sgv0Im[jroot])
                   + QSDdotc(qfncIm[jroot], sgv0Im[iroot]);
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

  // 0 quantum number
  Qshapes<> qz(1, Quantum::zero());
  Dshapes   dz(qz.size(), 1);

  TVector<Qshapes<>, 3>qshape = make_array( qz, qz, qz);
  TVector<Dshapes,   3>dshape = make_array( dz, dz, dz);

  for(int i = 0; i < L; ++i) {
    for(int iroot = 0; iroot < sites[i].nstore; ++iroot) {
      sites[i].lopr[iroot].resize(Quantum::zero(), qshape, dshape);
      sites[i].lopr[iroot] = 0.0;
      sites[i].ropr[iroot].resize(Quantum::zero(), qshape, dshape);
      sites[i].ropr[iroot] = 0.0;
    }
    sites[i].save(i);
  }
}

int prototype::lrt_sweep(bool forward, MpStorages& sites0, MpStorages& sitesRe, MpStorages& sitesIm,
                         DArray<2>& a_subspace, DArray<2>& b_subspace, DArray<2>& s_subspace, DArray<2>& d_subspace,
                         vector<double>& eigv, vector<double>& rnorm, vector<double>& bnorm, int nroots, int mroots, int kroots, const CALC_TYPE& calc_type)
{
  DArray<2> alphaRe;
  DArray<2> alphaIm;

  if(mroots > 1) {
    printing_matrix_elements(a_subspace, b_subspace, s_subspace, d_subspace);

//  DEBUG:
//  ComputeEigenvaluesTest(a_subspace, b_subspace, s_subspace, d_subspace);
//  DEBUG:
    DArray<1> eigvtmp;
    int reduced_mroots = 1 + ComputeEigenvalues(a_subspace, b_subspace, s_subspace, d_subspace, eigvtmp, alphaRe, alphaIm);

    if(reduced_mroots < mroots)
      mroots = reduced_mroots;

    cout << "\tEigenvalues" << endl;
    for(int i = 1; i < nroots; ++i) {
      eigv[i] = eigvtmp(i-1);
      cout << "\t\tState [ " << setw(2) << i << " ]: " << setw(16) << fixed << setprecision(12) << eigv[i] << endl;
    }
  }

  int lroots = mroots + nroots - kroots;

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

//RotateStorage(alphaRe, sitesRe[0].wfnc, mroots);
//RotateStorage(alphaIm, sitesIm[0].wfnc, mroots);
  RotateStorage(alphaRe, alphaIm, sitesRe[0].wfnc, sitesIm[0].wfnc, mroots);

  bool boundary = true;

  for(int i = 0; i < L-1; ++i) {
    cout << "\t\t\tSTARTING BLOCK ITERATION FOR SITE " << i << endl;
    cout << "\t====================================================================================================" << endl;
    sites0 [i+1].load(i+1);
    sitesRe[i+1].load(i+1);
    sitesIm[i+1].load(i+1);

    QSDArray<3> diag(sites0[i].wfnc[0].q(), sites0[i].wfnc[0].qshape());
    ComputeDiagonal(sites0[i].mpo, sites0[i].lopr[0], sites0[i].ropr[0], diag);
    for(QSDArray<3>::iterator it = diag.begin(); it != diag.end(); ++it) {
      for(DArray<3>::iterator id = it->second->begin(); id != it->second->end(); ++id) {
        *id -= eigv[0];
      }
    }

//  RotateStorage(alphaRe, sitesRe[i].ropr, mroots);
//  RotateStorage(alphaIm, sitesIm[i].ropr, mroots);
    RotateStorage(alphaRe, alphaIm, sitesRe[i].ropr, sitesIm[i].ropr, mroots);

    cout << "\tSolving correction equation @ site " << i << endl;
    SolveCorrectionEquation(1, sites0[i], sitesRe[i], sitesIm[i], diag, eigv, rnorm, nroots, mroots, kroots, boundary, calc_type);

//  RotateStorage(alphaRe, sitesRe[i+1].rmps, mroots);
//  RotateStorage(alphaIm, sitesIm[i+1].rmps, mroots);
    RotateStorage(alphaRe, alphaIm, sitesRe[i+1].rmps, sitesIm[i+1].rmps, mroots);

    for(int iroot = 1; iroot < mroots; ++iroot) {
      sitesRe[i+1].wfnc[iroot].clear();
      ComputeGuess(1, sites0[i].lmps[0], sites0[i].wfnc[0], sitesRe[i+1].rmps[iroot], sitesRe[i+1].wfnc[iroot]);
      ComputeGuess(1, sites0[i].lmps[0], sitesRe[i].wfnc[iroot], sites0[i+1].rmps[0], sitesRe[i+1].wfnc[iroot]);
      Orthogonalize(sites0[i+1].wfnc[0], sitesRe[i+1].wfnc[iroot]);

      sitesIm[i+1].wfnc[iroot].clear();
      ComputeGuess(1, sites0[i].lmps[0], sites0[i].wfnc[0], sitesIm[i+1].rmps[iroot], sitesIm[i+1].wfnc[iroot]);
      ComputeGuess(1, sites0[i].lmps[0], sitesIm[i].wfnc[iroot], sites0[i+1].rmps[0], sitesIm[i+1].wfnc[iroot]);
      Orthogonalize(sites0[i+1].wfnc[0], sitesIm[i+1].wfnc[iroot]);
    }
//  for(int iroot = mroots; iroot < lroots; ++iroot) {
//    sitesRe[i+1].wfnc[iroot].clear();
//    ComputeGuess(1, sites0[i].lmps[0], sitesRe[i].wfnc[iroot], sites0[i+1].rmps[0], sitesRe[i+1].wfnc[iroot]);
//    Orthogonalize(sites0[i+1].wfnc[0], sitesRe[i+1].wfnc[iroot]);

//    sitesIm[i+1].wfnc[iroot].clear();
//    ComputeGuess(1, sites0[i].lmps[0], sitesIm[i].wfnc[iroot], sites0[i+1].rmps[0], sitesIm[i+1].wfnc[iroot]);
//    Orthogonalize(sites0[i+1].wfnc[0], sitesIm[i+1].wfnc[iroot]);
//  }

//-- to be cut
//  cout << "\tComputing matrix elements @ site " << i << endl;
//  ComputeMatrixElements(1, sites0[i], sitesRe[i], sitesIm[i], eigv, h_subspace, s_subspace, bnorm, nroots, mroots, kroots, boundary, false);
//  printing_matrix_elements(h_subspace, s_subspace);
//--

    boundary = false;

    QSDArray<2> gauge_0;
    QSDgemm(ConjTrans, NoTrans, 1.0, sites0[i].lmps[0], sites0[i].wfnc[0], 1.0, gauge_0);
    QSDArray<2> gauge_i;
    ComputeInverseGauge(gauge_0, gauge_i);

//QSDArray<2> check0i;
//QSDgemm(NoTrans, NoTrans, 1.0, gauge_0, gauge_i, 1.0, check0i);
//cout << "debug: gauge-trans 0-i check = " << check0i << endl;
//QSDArray<2> checki0;
//QSDgemm(NoTrans, NoTrans, 1.0, gauge_i, gauge_0, 1.0, checki0);
//cout << "debug: gauge-trans i-0 check = " << checki0 << endl;
//  checking_gauge_condition(sites0[i], sitesRe[i], sites0[i+1], sitesRe[i+1], 1);

    for(int iroot = 1; iroot < lroots; ++iroot) {
      QSDArray<3> tmpRe(sitesRe[i].wfnc[iroot]);
      Orthogonalize(1, sites0[i].lmps[0], tmpRe);
      sitesRe[i].lmps[iroot].clear();
      QSDgemm(NoTrans, NoTrans, 1.0, tmpRe, gauge_i, 1.0, sitesRe[i].lmps[iroot]);

      sitesRe[i+1].lopr[iroot].clear();
      Renormalize(1, sites0[i].mpo, sitesRe[i].lopr[iroot], sites0[i].lmps[0], sites0[i].lmps[0], sitesRe[i+1].lopr[iroot]);
      Renormalize(1, sites0[i].mpo, sites0[i].lopr[0], sites0[i].lmps[0], sitesRe[i].lmps[iroot], sitesRe[i+1].lopr[iroot]);

      QSDArray<3> tmpIm(sitesIm[i].wfnc[iroot]);
      Orthogonalize(1, sites0[i].lmps[0], tmpIm);
      sitesIm[i].lmps[iroot].clear();
      QSDgemm(NoTrans, NoTrans, 1.0, tmpIm, gauge_i, 1.0, sitesIm[i].lmps[iroot]);

      sitesIm[i+1].lopr[iroot].clear();
      Renormalize(1, sites0[i].mpo, sitesIm[i].lopr[iroot], sites0[i].lmps[0], sites0[i].lmps[0], sitesIm[i+1].lopr[iroot]);
      Renormalize(1, sites0[i].mpo, sites0[i].lopr[0], sites0[i].lmps[0], sitesIm[i].lmps[iroot], sitesIm[i+1].lopr[iroot]);
    }

    sites0 [i].save(i);
    sitesRe[i].save(i);
    sitesIm[i].save(i);
    cout << "\t====================================================================================================" << endl;
  }

  {
    QSDArray<3> diag(sites0[L-1].wfnc[0].q(), sites0[L-1].wfnc[0].qshape());
    ComputeDiagonal(sites0[L-1].mpo, sites0[L-1].lopr[0], sites0[L-1].ropr[0], diag);
    for(QSDArray<3>::iterator it = diag.begin(); it != diag.end(); ++it) {
      for(DArray<3>::iterator id = it->second->begin(); id != it->second->end(); ++id) {
        *id -= eigv[0];
      }
    }

//  RotateStorage(alphaRe, sitesRe[L-1].ropr, mroots);
//  RotateStorage(alphaIm, sitesIm[L-1].ropr, mroots);
    RotateStorage(alphaRe, alphaIm, sitesRe[L-1].ropr, sitesIm[L-1].ropr, mroots);

    cout << "\tSolving correction equation @ site " << L-1 << endl;
    SolveCorrectionEquation(1, sites0[L-1], sitesRe[L-1], sitesIm[L-1], diag, eigv, rnorm, nroots, mroots, kroots, boundary, calc_type);
//-- to be cut
//  cout << "\tComputing matrix elements @ site " << L-1 << endl;
//  ComputeMatrixElements(1, sites0[L-1], sitesRe[L-1], sitesIm[L-1], eigv, h_subspace, s_subspace, bnorm, nroots, mroots, kroots, true, true);
//  printing_matrix_elements(h_subspace, s_subspace);
//--

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

//-- to be cut
//RotateStorage(alpha, sitesRe[L-1].wfnc, mroots);
//RotateStorage(alpha, sitesIm[L-1].wfnc, mroots);
//--

  bool boundary = true;

  for(int i = L-1; i > 0; --i) {
    cout << "\t\t\tSTARTING BLOCK ITERATION FOR SITE " << i << endl;
    cout << "\t====================================================================================================" << endl;
    sites0 [i-1].load(i-1);
    sitesRe[i-1].load(i-1);
    sitesIm[i-1].load(i-1);

//-- to be cut
//  QSDArray<3> diag(sites0[i].wfnc[0].q(), sites0[i].wfnc[0].qshape());
//  ComputeDiagonal(sites0[i].mpo, sites0[i].lopr[0], sites0[i].ropr[0], diag);
//  for(QSDArray<3>::iterator it = diag.begin(); it != diag.end(); ++it) {
//    for(DArray<3>::iterator id = it->second->begin(); id != it->second->end(); ++id) {
//      *id -= eigv[0];
//    }
//  }

//  RotateStorage(alpha, sitesRe[i].lopr, mroots);
//  RotateStorage(alpha, sitesIm[i].lopr, mroots);

//  cout << "\tSolving correction equation @ site " << i << endl;
//  SolveCorrectionEquation(sites0[i], sitesRe[i], sitesIm[i], diag, eigv, rnorm, nroots, mroots, kroots, boundary);

//  RotateStorage(alpha, sitesRe[i-1].lmps, mroots);
//  RotateStorage(alpha, sitesIm[i-1].lmps, mroots);
//--

//  for(int iroot = 1; iroot < mroots; ++iroot) {
    for(int iroot = 1; iroot < lroots; ++iroot) {
      sitesRe[i-1].wfnc[iroot].clear();
      ComputeGuess(0, sites0[i].rmps[0], sites0[i].wfnc[0], sitesRe[i-1].lmps[iroot], sitesRe[i-1].wfnc[iroot]);
      ComputeGuess(0, sites0[i].rmps[0], sitesRe[i].wfnc[iroot], sites0[i-1].lmps[0], sitesRe[i-1].wfnc[iroot]);
      Orthogonalize(sites0[i-1].wfnc[0], sitesRe[i-1].wfnc[iroot]);

      sitesIm[i-1].wfnc[iroot].clear();
      ComputeGuess(0, sites0[i].rmps[0], sites0[i].wfnc[0], sitesIm[i-1].lmps[iroot], sitesIm[i-1].wfnc[iroot]);
      ComputeGuess(0, sites0[i].rmps[0], sitesIm[i].wfnc[iroot], sites0[i-1].lmps[0], sitesIm[i-1].wfnc[iroot]);
      Orthogonalize(sites0[i-1].wfnc[0], sitesIm[i-1].wfnc[iroot]);
    }

    cout << "\tComputing matrix elements @ site " << i << endl;
    ComputeMatrixElements(0, sites0[i], sitesRe[i], sitesIm[i], eigv, a_subspace, b_subspace, s_subspace, d_subspace, bnorm, nroots, mroots, kroots, boundary, false, calc_type);
//  printing_matrix_elements(a_subspace, b_subspace, s_subspace, d_subspace);

    boundary = false;

    QSDArray<2> gauge_0;
    QSDgemm(NoTrans, ConjTrans, 1.0, sites0[i].wfnc[0], sites0[i].rmps[0], 1.0, gauge_0);
    QSDArray<2> gauge_i;
    ComputeInverseGauge(gauge_0, gauge_i);

    for(int iroot = 1; iroot < lroots; ++iroot) {
      QSDArray<3> tmpRe(sitesRe[i].wfnc[iroot]);
      Orthogonalize(0, sites0[i].rmps[0], tmpRe);
      sitesRe[i].rmps[iroot].clear();
      QSDgemm(NoTrans, NoTrans, 1.0, gauge_i, tmpRe, 1.0, sitesRe[i].rmps[iroot]);

      sitesRe[i-1].ropr[iroot].clear();
      Renormalize(0, sites0[i].mpo, sitesRe[i].ropr[iroot], sites0[i].rmps[0], sites0[i].rmps[0], sitesRe[i-1].ropr[iroot]);
      Renormalize(0, sites0[i].mpo, sites0[i].ropr[0], sites0[i].rmps[0], sitesRe[i].rmps[iroot], sitesRe[i-1].ropr[iroot]);

      QSDArray<3> tmpIm(sitesIm[i].wfnc[iroot]);
      Orthogonalize(0, sites0[i].rmps[0], tmpIm);
      sitesIm[i].rmps[iroot].clear();
      QSDgemm(NoTrans, NoTrans, 1.0, gauge_i, tmpIm, 1.0, sitesIm[i].rmps[iroot]);

      sitesIm[i-1].ropr[iroot].clear();
      Renormalize(0, sites0[i].mpo, sitesIm[i].ropr[iroot], sites0[i].rmps[0], sites0[i].rmps[0], sitesIm[i-1].ropr[iroot]);
      Renormalize(0, sites0[i].mpo, sites0[i].ropr[0], sites0[i].rmps[0], sitesIm[i].rmps[iroot], sitesIm[i-1].ropr[iroot]);
    }

    sites0 [i].save(i);
    sitesRe[i].save(i);
    sitesIm[i].save(i);
    cout << "\t====================================================================================================" << endl;
  }

  {
//-- to be cut
//  QSDArray<3> diag(sites0[0].wfnc[0].q(), sites0[0].wfnc[0].qshape());
//  ComputeDiagonal(sites0[0].mpo, sites0[0].lopr[0], sites0[0].ropr[0], diag);
//  for(QSDArray<3>::iterator it = diag.begin(); it != diag.end(); ++it) {
//    for(DArray<3>::iterator id = it->second->begin(); id != it->second->end(); ++id) {
//      *id -= eigv[0];
//    }
//  }

//  RotateStorage(alpha, sitesRe[0].lopr, mroots);
//  RotateStorage(alpha, sitesIm[0].lopr, mroots);

//  cout << "\tSolving correction equation @ site " << 0 << endl;
//  SolveCorrectionEquation(sites0[0], sitesRe[0], sitesIm[0], diag, eigv, rnorm, nroots, mroots, kroots, true);
//--
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

  int mxritz = max(3*(nroots-1), 20);
  int nstore = (nroots == 1) ? 1 : (nroots+mxritz);

  MpStorages sitesRe(L, MpSite("state-Re", nstore, input.prefix));
  StorageInitialize(sitesRe);

  MpStorages sitesIm(L, MpSite("state-Im", nstore, input.prefix));
  StorageInitialize(sitesIm);

  DArray<2> a_subspace;
  DArray<2> b_subspace;
  DArray<2> s_subspace;
  DArray<2> d_subspace;

  vector<double> eigvs(nroots, 0.0);
  vector<double> eigvc(nroots, 0.0);
  vector<double> rnorm;
  vector<double> bnorm;

  eigvs[0] = input.energy;
  eigvc[0] = input.energy;
  bool forward = 1;
  int nconv_roots = 1;
  int mroots = 1;

  //
  // linear-response
  //

  // TDA

  for(int iter = 0; iter < 100; ++iter) {
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

    mroots = lrt_sweep(1, sites0, sitesRe, sitesIm, a_subspace, b_subspace, s_subspace, d_subspace, eigvs, rnorm, bnorm, nroots, mroots, nconv_roots, TDA);
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
        if(fabs(rnorm[nconv_roots]) >= input.tolerance*10) break;
    }

    if(nconv_roots == nroots) break;

//  mroots += nroots - nconv_roots;
    if(mroots > mxritz) mroots = nroots;

// FOR TEST (don't perform correction)
//  mroots = nroots;
//  nconv_roots = nroots;
  }

  // RPA

  mroots = 1;
  nconv_roots = 1;

  for(int iter = 0; iter < 100; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tDAVIDSON SWEEP ITERATION [ " << setw(4) << iter << " ] FOR RPA"   << endl;
    cout << "\t====================================================================================================" << endl;

    mroots = lrt_sweep(1, sites0, sitesRe, sitesIm, a_subspace, b_subspace, s_subspace, d_subspace, eigvc, rnorm, bnorm, nroots, mroots, nconv_roots, RPA);

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
        if(fabs(rnorm[nconv_roots]) >= input.tolerance*10) break;
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
    ecorr1 -= eigvs[iroot] * bnorm[iroot];
    ecorr2 += 0.5 * (eigvc[iroot] - eigvs[iroot]);
  }

  cout.precision(16);
  cout << "\t====================================================================================================" << endl;
  cout << "\tDEBUG: Correlation Energy (Scheme 1) = " << setw(24) << fixed << ecorr1 << " ( used ) " << endl;
  cout << "\tDEBUG: Correlation Energy (Scheme 2) = " << setw(24) << fixed << ecorr2 << endl;
  double ecorr = ecorr1;
  cout << "\t====================================================================================================" << endl;
  cout << "\t\tPRINTING RPA ENERGY" << endl;
  cout.precision(16);
  cout << "\t\t\tEnergy[ 0] = " << setw(24) << fixed << eigvs[0] + ecorr << " ( E(corr) = " << setw(24) << fixed << ecorr << " ) " << endl;
  for(int i = 1; i < nroots; ++i) {
    cout.precision(16);
    cout << "\t\t\tEnergy[" << setw(2) << i << "] = " << setw(24) << fixed << eigvs[i] + eigvs[0] + ecorr << endl;
  }

  return eigvs;
}
