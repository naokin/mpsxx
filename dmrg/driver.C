#include <iostream>
#include <iomanip>
using namespace std;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; };

#include <btas/SPARSE/SDdiagonal.h>
#include <btas/QSPARSE/QSDArray.h>
#include "btas_template_specialize.h"

#include "driver.h"
using btas::shape;
using btas::NoTrans;
using btas::Trans;
using btas::ConjTrans;

void prototype::ComputeGuess
(bool forward, const btas::QSDArray<3>& mps0,
               const btas::QSDArray<3>& wfn0,
               const btas::QSDArray<3>& mps1,
                     btas::QSDArray<3>& wfn1)
{
  using btas::NoTrans;
  using btas::ConjTrans;

  if(forward) {
    btas::QSDArray<2> lres;
    btas::QSDgemm(ConjTrans, NoTrans, 1.0, mps0, wfn0, 1.0, lres);
    btas::QSDgemm(  NoTrans, NoTrans, 1.0, lres, mps1, 1.0, wfn1);
  }
  else {
    btas::QSDArray<2> rres;
    btas::QSDgemm(NoTrans, ConjTrans, 1.0, wfn0, mps0, 1.0, rres);
    btas::QSDgemm(NoTrans,   NoTrans, 1.0, mps1, rres, 1.0, wfn1);
  }
}

void prototype::Canonicalize
(bool forward, const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& mps0, int M)
{
  if(forward) {
    btas::SDArray<1> s;
    btas::QSDArray<2> v;
    btas::QSDgesvd(btas::LeftArrow,  wfn0, s, mps0, v, M);
  }
  else {
    btas::SDArray<1> s;
    btas::QSDArray<2> u;
    btas::QSDgesvd(btas::RightArrow, wfn0, s, u, mps0, M);
  }
}

void prototype::Canonicalize
(bool forward, const btas::QSDArray<4>& wfnx,
                     btas::QSDArray<3>& mps0,
                     btas::QSDArray<3>& wfn1, int M)
{
  if(forward) {
    btas::SDArray <1> s;
    btas::QSDgesvd(btas::LeftArrow,  wfnx, s, mps0, wfn1, M);
    btas::SDdidm(s, wfn1);
  }
  else {
    btas::SDArray <1> s;
    btas::QSDgesvd(btas::RightArrow, wfnx, s, wfn1, mps0, M);
    btas::SDdimd(wfn1, s);
  }
}

void prototype::Renormalize
(bool forward, const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& opr0,
               const btas::QSDArray<3>& bra0,
               const btas::QSDArray<3>& ket0,
                     btas::QSDArray<3>& opr1)
{
  if(forward) {
    btas::QSDArray<4> scr1;
    btas::QSDcontract(1.0, opr0, shape(0), bra0.conjugate(), shape(0), 1.0, scr1);
    btas::QSDArray<4> scr2;
    btas::QSDcontract(1.0, scr1, shape(0, 2), mpo0, shape(0, 1), 1.0, scr2);
    btas::QSDcontract(1.0, scr2, shape(0, 2), ket0, shape(0, 1), 1.0, opr1);
  }
  else {
    btas::QSDArray<4> scr1;
    btas::QSDcontract(1.0, bra0.conjugate(), shape(2), opr0, shape(0), 1.0, scr1);
    btas::QSDArray<4> scr2;
    btas::QSDcontract(1.0, scr1, shape(1, 2), mpo0, shape(1, 3), 1.0, scr2);
    btas::QSDcontract(1.0, scr2, shape(3, 1), ket0, shape(1, 2), 1.0, opr1);
  }
}

void prototype::ComputeInverseGauge
(              const btas::QSDArray<2>& gauge_0,
                     btas::QSDArray<2>& gauge_i)
{
  btas::SDArray<1> s;
  btas::QSDArray<2> u;
  btas::QSDArray<2> v;
  btas::QSDgesvd(btas::LeftArrow, gauge_0, s, u, v, 0);
  for(btas::SDArray<1>::iterator it = s.begin(); it != s.end(); ++it) {
    for(btas::DArray<1>::iterator id = it->second->begin(); id != it->second->end(); ++id) {
      if(fabs(*id) >= 1.0e-20)
        *id = 1.0/(*id);
      else
        *id = 0.0;
    }
  }
  btas::SDdimd(u, s);
  gauge_i.clear();
  btas::QSDgemm(btas::ConjTrans, btas::ConjTrans, 1.0, v, u, 1.0, gauge_i);
}

void prototype::Normalize(btas::QSDArray<3>& wfn0)
{
  double norm = btas::QSDdotc(wfn0, wfn0);
  if(norm >= 1.0e-20) {
    btas::QSDscal(1.0/sqrt(norm), wfn0);
  }
}

void prototype::Orthogonalize(const btas::QSDArray<3>& wfn0, btas::QSDArray<3>& wfn1)
{
  double overlap = btas::QSDdotc(wfn0, wfn1);
  QSDaxpy(-overlap, wfn0, wfn1);
}

void prototype::Orthogonalize(bool forward, const btas::QSDArray<3>& mps0, btas::QSDArray<3>& wfn1)
{
  if(forward) {
    btas::QSDArray<2> proj;
    btas::QSDgemm(ConjTrans, NoTrans, 1.0, mps0, wfn1, 1.0, proj);
    btas::QSDgemm(  NoTrans, NoTrans,-1.0, mps0, proj, 1.0, wfn1);
  }
  else {
    btas::QSDArray<2> proj;
    btas::QSDgemm(NoTrans, ConjTrans, 1.0, wfn1, mps0, 1.0, proj);
    btas::QSDgemm(NoTrans,   NoTrans,-1.0, proj, mps0, 1.0, wfn1);
  }
}

void prototype::ComputeDiagonal
(              const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
                     btas::QSDArray<3>& diag)
{
  btas::SDArray<3> mpo0_diag;
  btas::SDArray<2> lopr_diag;
  btas::SDArray<2> ropr_diag;

  btas::SDdiagonal(mpo0, shape(1, 2), mpo0_diag);
  btas::SDdiagonal(lopr, shape(0, 2), lopr_diag);
  btas::SDdiagonal(ropr, shape(0, 2), ropr_diag);

  btas::SDArray<3> scr1;
  btas::SDcontract(1.0, lopr_diag, shape(1), mpo0_diag, shape(0), 1.0, scr1);
  btas::SDArray<3> scr2;
  btas::SDcontract(1.0, scr1,      shape(2), ropr_diag, shape(1), 1.0, scr2);
  btas::SDcopy(scr2, diag, 1); // preserve quantum number of diag
}

void prototype::ComputeDiagonal
(              const btas::QSDArray<4>& lmpo,
               const btas::QSDArray<4>& rmpo,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
                     btas::QSDArray<4>& diag)
{
  btas::SDArray<3> lmpo_diag;
  btas::SDArray<3> rmpo_diag;
  btas::SDArray<2> lopr_diag;
  btas::SDArray<2> ropr_diag;

  btas::SDdiagonal(lmpo, shape(1, 2), lmpo_diag);
  btas::SDdiagonal(rmpo, shape(1, 2), rmpo_diag);
  btas::SDdiagonal(lopr, shape(0, 2), lopr_diag);
  btas::SDdiagonal(ropr, shape(0, 2), ropr_diag);

  btas::SDArray<3> scr1;
  btas::SDcontract(1.0, lopr_diag, shape(1), lmpo_diag, shape(0), 1.0, scr1);
  btas::SDArray<4> scr2;
  btas::SDcontract(1.0, scr1,      shape(2), rmpo_diag, shape(0), 1.0, scr2);
  btas::SDArray<4> scr3;
  btas::SDcontract(1.0, scr2,      shape(3), ropr_diag, shape(1), 1.0, scr3);
  btas::SDcopy(scr3, diag, 1); // preserve quantum number of diag
}

void prototype::ComputeSigmaVector // compute: H | psi >
(              const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
               const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& sgv0)
{
  btas::QSDArray<4> scr1;
  btas::QSDcontract(1.0, lopr, shape(2), wfn0, shape(0), 1.0, scr1);
  btas::QSDArray<4> scr2;
  btas::QSDcontract(1.0, scr1, shape(1, 2), mpo0, shape(0, 2), 1.0, scr2);
  btas::QSDcontract(1.0, scr2, shape(3, 1), ropr, shape(1, 2), 1.0, sgv0);
}

void prototype::ComputeSigmaVectorConj // compute: < psi | H  ... but returned as ket state | sgv >
(              const btas::QSDArray<4>& mpo0,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
               const btas::QSDArray<3>& wfn0,
                     btas::QSDArray<3>& sgv0)
{
  sgv0.conjugate_self();
  btas::QSDArray<4> scr1;
  btas::QSDcontract(1.0, lopr, shape(0), wfn0.conjugate(), shape(0), 1.0, scr1);
  btas::QSDArray<4> scr2;
  btas::QSDcontract(1.0, scr1, shape(0, 2), mpo0, shape(0, 1), 1.0, scr2);
  btas::QSDcontract(1.0, scr2, shape(1, 3), ropr, shape(0, 1), 1.0, sgv0);
  sgv0.conjugate_self();
}

void prototype::ComputeSigmaVector
(              const btas::QSDArray<4>& lmpo,
               const btas::QSDArray<4>& rmpo,
               const btas::QSDArray<3>& lopr,
               const btas::QSDArray<3>& ropr,
               const btas::QSDArray<4>& wfn0,
                     btas::QSDArray<4>& sgv0)
{
  btas::QSDArray<5> scr1;
  btas::QSDcontract(1.0, lopr, shape(2), wfn0, shape(0), 1.0, scr1);
  btas::QSDArray<5> scr2;
  btas::QSDcontract(1.0, scr1, shape(1, 2), lmpo, shape(0, 2), 1.0, scr2);
  btas::QSDArray<5> scr3;
  btas::QSDcontract(1.0, scr2, shape(4, 1), rmpo, shape(0, 2), 1.0, scr3);
  btas::QSDcontract(1.0, scr3, shape(4, 1), ropr, shape(1, 2), 1.0, sgv0);
}

int prototype::ComputeOrthogonalTransform
(              const btas::DArray<2>& s_matrix,
                     btas::DArray<2>& u_matrix)
{
  int nrows = s_matrix.shape(0);
  // diagonalize s_subspace
  btas::DArray<1> d;
  btas::DArray<2> p_matrix;
  btas::Dsyev(s_matrix, d, p_matrix);

//cout << "\tprinting eigenvalues of overlap matrix" << endl;
//for(int i = 0; i < d.size(); ++i)
//  cout << "\t\tState [ " << setw(2) << i + 1 << " ]: " << setw(16) << fixed << setprecision(12) << d(i) << endl;

  // check positive eigenvalues;
  int nimag = 0;
  for(; nimag < d.size(); ++nimag) if(d(nimag) >= 1.0e-20) break;

  int ncols = nrows - nimag;

  btas::DArray<1> di(ncols);
  u_matrix.resize(nrows, ncols);

  int icol = 0;
  for(int i = nimag; i < p_matrix.shape(1); ++i, ++icol) {
    di(icol) = 1.0/sqrt(d(i));
    for(int j = 0; j < p_matrix.shape(0); ++j) {
      u_matrix(j, icol) = p_matrix(j, i);
    }
  }
  btas::Ddimd(u_matrix, di);

  return ncols;
}

void prototype::TransformedMatrix
(              const btas::DArray<2>& h_matrix,
               const btas::DArray<2>& u_matrix,
                     btas::DArray<2>& h_transm)
{
  btas::DArray<2> h_tmpmat;
  btas::Dgemm(NoTrans, NoTrans, 1.0, h_matrix, u_matrix, 1.0, h_tmpmat);
  btas::Dgemm(  Trans, NoTrans, 1.0, u_matrix, h_tmpmat, 1.0, h_transm);
}

int prototype::ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& s_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alpha)
{
  btas::DArray<2> u_matrix;
  int ncols = ComputeOrthogonalTransform(s_subspace, u_matrix);

  btas::DArray<2> a_trans;
  TransformedMatrix(a_subspace, u_matrix, a_trans);

  btas::DArray<2> t_matrix;
  btas::Dsyev(a_trans, eigvs, t_matrix);
  alpha.clear();
//btas::Dgemm(NoTrans, Trans, 1.0, u_matrix, t_matrix, 1.0, alpha);
  btas::Dgemm(NoTrans, NoTrans, 1.0, u_matrix, t_matrix, 1.0, alpha);

  return ncols;
}

/*
int prototype::ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm)
{
  int nrows = a_subspace.shape(0);
  int ncols = a_subspace.shape(1);

  btas::DArray<2> s_matrix(2*nrows, 2*ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      s_matrix(i,       j      ) =  s_subspace(i, j);
      s_matrix(i,       j+ncols) =  d_subspace(i, j);
      s_matrix(i+nrows, j      ) = -d_subspace(i, j);
      s_matrix(i+nrows, j+ncols) = -s_subspace(i, j);
    }
  }
  btas::DArray<2> h_matrix(2*nrows, 2*ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      h_matrix(i,       j      ) =  a_subspace(i, j);
      h_matrix(i+nrows, j+ncols) =  a_subspace(i, j);
      h_matrix(i,       j+ncols) =  b_subspace(i, j);
      h_matrix(i+nrows, j      ) =  b_subspace(i, j);
    }
  }

  btas::DArray<2> u_matrix;
  ncols = ComputeOrthogonalTransform(s_matrix, u_matrix);

  btas::DArray<2> h_reduce;
  TransformedMatrix(h_matrix, u_matrix, h_reduce);

  btas::DArray<1> d;
  btas::DArray<2> p;
  btas::Dsyev(h_reduce, d, p);

  // check positive eigenvalues;
  int nimag = 0;
  cout.precision(8);
  for(; nimag < d.size(); ++nimag) {
    if(d(nimag) >= 1.0e-12) break;
    cout << "\tImaginary Eigenvalue " << setw(12) << fixed << d(nimag) << " is found @ state [" << setw(2) << nimag << "] " << endl;
  }

  btas::DArray<2> x_matrix(ncols-nimag, ncols);
  eigvs.resize(ncols-nimag);

  int icol = 0;
  for(int i = nimag; i < ncols; ++i, ++icol) {
    eigvs(icol) = d(i);
    for(int j = 0; j < ncols; ++j) {
      x_matrix(icol, j) = p(i, j);
    }
  }
  ncols -= nimag;

  btas::DArray<2> alpha;
  btas::Dgemm(NoTrans, Trans, 1.0, u_matrix, x_matrix, 1.0, alpha);

  // DEBUG CHECK DIAGONAL
  btas::DArray<2> hreal;
  btas::DArray<2> hrealtmp;
  btas::Dgemm(NoTrans, NoTrans, 1.0, s_matrix, alpha, 1.0, hrealtmp);
  btas::Dgemm(  Trans, NoTrans, 1.0, alpha, hrealtmp, 1.0, hreal);
  cout << "DEBUG: H(real)" << endl;
  for(int i = 0; i < hreal.shape(0); ++i) {
    cout << "\t";
    for(int j = 0; j < hreal.shape(1); ++j) {
      cout << setw(12) << fixed << hreal(i, j);
    }
    cout << endl;
  }
  btas::DArray<2> alphaTr(2*nrows, ncols);

  alphaRe.resize(nrows, ncols);
  alphaIm.resize(nrows, ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      alphaRe(i, j) = alpha(i,       j);
      alphaIm(i, j) = alpha(i+nrows, j);
      // FOR DEBUG
      alphaTr(i, j) = alpha(i+nrows, j);
      alphaTr(i+nrows, j) = alpha(i, j);
    }
  }

  btas::DArray<2> himag;
  btas::DArray<2> himagtmp;
  btas::Dgemm(NoTrans, NoTrans, 1.0, s_matrix, alphaTr, 1.0, himagtmp);
  btas::Dgemm(  Trans, NoTrans, 1.0, alphaTr, himagtmp, 1.0, himag);
  cout << "DEBUG: H(imag)" << endl;
  for(int i = 0; i < himag.shape(0); ++i) {
    cout << "\t";
    for(int j = 0; j < himag.shape(1); ++j) {
      cout << setw(12) << fixed << himag(i, j);
    }
    cout << endl;
  }

  return ncols;
}
*/

int prototype::ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm)
{
  int nrows = a_subspace.shape(0);
  int ncols = a_subspace.shape(1);

  btas::DArray<2> s_matrix(2*nrows, 2*ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      s_matrix(i,       j      ) =  s_subspace(i, j);
      s_matrix(i,       j+ncols) =  d_subspace(i, j);
      s_matrix(i+nrows, j      ) = -d_subspace(i, j);
      s_matrix(i+nrows, j+ncols) = -s_subspace(i, j);
    }
  }
  btas::DArray<2> h_matrix(2*nrows, 2*ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      h_matrix(i,       j      ) =  a_subspace(i, j);
      h_matrix(i+nrows, j+ncols) =  a_subspace(i, j);
      h_matrix(i,       j+ncols) =  b_subspace(i, j);
      h_matrix(i+nrows, j      ) =  b_subspace(i, j);
    }
  }

  btas::DArray<2> u_matrix;
  ncols = ComputeOrthogonalTransform(h_matrix, u_matrix);

// DEBUG
//cout.precision(4);
//cout << "debug: u_matrix: [" << u_matrix.shape(0) << "," << u_matrix.shape(1) << "]" << endl;
//for(int i = 0; i < u_matrix.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < u_matrix.shape(1); ++j) {
//    cout << setw(8) << fixed << u_matrix(i, j);
//  }
//  cout << endl;
//}

  btas::DArray<2> s_reduce;
  TransformedMatrix(s_matrix, u_matrix, s_reduce);

  btas::DArray<1> d;
  btas::DArray<2> p;
  btas::Dsyev(s_reduce, d, p); // h_matrix as metric
//btas::Dsygv(s_matrix, h_matrix, d, p); // h_matrix as metric

  // check positive eigenvalues;
  int nmin = 0;
  for(; nmin < d.size(); ++nmin)
    if(d(nmin) >= 1.0e-16) break;
  int nmax = ncols-1;
  for(; nmax >= 0; --nmax)
    if(d(nmax) < 1.0e+8) break;
  ncols = nmax - nmin + 1;

  btas::DArray<2> x_matrix;
//btas::Dgemm(NoTrans, Trans, 1.0, u_matrix, p, 1.0, x_matrix);
  btas::Dgemm(NoTrans, NoTrans, 1.0, u_matrix, p, 1.0, x_matrix);

// DEBUG
//cout.precision(4);
//cout << "debug: x_matrix: [" << x_matrix.shape(0) << "," << x_matrix.shape(1) << "]" << endl;
//for(int i = 0; i < x_matrix.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < x_matrix.shape(1); ++j) {
//    cout << setw(8) << fixed << x_matrix(i, j);
//  }
//  cout << endl;
//}

  btas::DArray<2> alpha(2*nrows, ncols);
  eigvs.resize(ncols);

  int icol = 0;
  for(int i = nmax; i >= nmin; --i, ++icol) {
    eigvs(icol) = 1.0/d(i);
    for(int j = 0; j < 2*nrows; ++j) {
      alpha(j, icol) = x_matrix(j, i)/sqrt(d(i));
    }
  }

// DEBUG
//cout.precision(4);
//cout << "debug: alpha: [" << alpha.shape(0) << "," << alpha.shape(1) << "]" << endl;
//for(int i = 0; i < alpha.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < alpha.shape(1); ++j) {
//    cout << setw(8) << fixed << alpha(i, j);
//  }
//  cout << endl;
//}

//// DEBUG CHECK DIAGONAL
//btas::DArray<2> h_diag;
//btas::DArray<2> h_diagtmp;
//btas::Dgemm(NoTrans, NoTrans, 1.0, h_matrix, alpha, 1.0, h_diagtmp);
//btas::Dgemm(  Trans, NoTrans, 1.0, alpha, h_diagtmp, 1.0, h_diag);
//cout << "DEBUG: H(diagonalized)" << endl;
//for(int i = 0; i < h_diag.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < h_diag.shape(1); ++j) {
//    cout << setw(8) << fixed << h_diag(i, j);
//  }
//  cout << endl;
//}
//btas::DArray<2> s_diag;
//btas::DArray<2> s_diagtmp;
//btas::Dgemm(NoTrans, NoTrans, 1.0, s_matrix, alpha, 1.0, s_diagtmp);
//btas::Dgemm(  Trans, NoTrans, 1.0, alpha, s_diagtmp, 1.0, s_diag);
//cout << "DEBUG: S(diagonalized)" << endl;
//for(int i = 0; i < s_diag.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < s_diag.shape(1); ++j) {
//    cout << setw(8) << fixed << s_diag(i, j);
//  }
//  cout << endl;
//}
//// DEBUG

  alphaRe.resize(nrows, ncols);
  alphaIm.resize(nrows, ncols);
  for(int i = 0; i < nrows; ++i) {
    for(int j = 0; j < ncols; ++j) {
      alphaRe(i, j) = alpha(i,       j);
      alphaIm(i, j) = alpha(i+nrows, j);
    }
  }

// DEBUG
//cout.precision(4);
//cout << "debug: alphaRe: [" << alphaRe.shape(0) << "," << alphaRe.shape(1) << "]" << endl;
//for(int i = 0; i < alphaRe.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < alphaRe.shape(1); ++j) {
//    cout << setw(8) << fixed << alphaRe(i, j);
//  }
//  cout << endl;
//}
//cout << "debug: alphaIm: [" << alphaIm.shape(0) << "," << alphaIm.shape(1) << "]" << endl;
//for(int i = 0; i < alphaIm.shape(0); ++i) {
//  cout << "\t";
//  for(int j = 0; j < alphaIm.shape(1); ++j) {
//    cout << setw(8) << fixed << alphaIm(i, j);
//  }
//  cout << endl;
//}

  return ncols;
}

void prototype::ComputeEigenvaluesTest
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace)
{
  int nrows = a_subspace.shape(0);
  int ncols = a_subspace.shape(1);
  btas::DArray<2> p_matrix;
  btas::DArray<2> m_matrix;
  btas::Dcopy(s_subspace, p_matrix);
  btas::Dcopy(s_subspace, m_matrix);
  btas::Daxpy( 1.0, d_subspace, p_matrix);
  btas::Daxpy(-1.0, d_subspace, m_matrix);

  btas::DArray<2> pxm;
  btas::Dgemm(NoTrans, NoTrans, 1.0, p_matrix, m_matrix, 1.0, pxm);

  btas::DArray<1> dpxm;
  btas::DArray<2> upxm;
  btas::Dsyev(pxm, dpxm, upxm);

  btas::DArray<2> mxp;
  btas::Dgemm(NoTrans, NoTrans, 1.0, m_matrix, p_matrix, 1.0, mxp);

  btas::DArray<1> dmxp;
  btas::DArray<2> umxp;
  btas::Dsyev(mxp, dmxp, umxp);

  cout.precision(8);
  cout << "DEBUG: (S + D)(S - D)" << endl;
  for(int i = 0; i < pxm.shape(0); ++i) {
    cout << "\t";
    for(int j = 0; j < pxm.shape(1); ++j) {
      cout << setw(12) << fixed << pxm(i, j);
    }
    cout << endl;
  }
  cout << "\tEigenvalues" << endl;
  for(int i = 0; i < dpxm.size(); ++i) {
    cout << "\t" << setw(12) << fixed << dpxm(i) << endl;
  }

  cout << "DEBUG: (S - D)(S + D)" << endl;
  for(int i = 0; i < mxp.shape(0); ++i) {
    cout << "\t";
    for(int j = 0; j < mxp.shape(1); ++j) {
      cout << setw(12) << fixed << mxp(i, j);
    }
    cout << endl;
  }
  cout << "\tEigenvalues" << endl;
  for(int i = 0; i < dmxp.size(); ++i) {
    cout << "\t" << setw(12) << fixed << dmxp(i) << endl;
  }
}

void prototype::AnalyzeTransferOperator
(              const btas::QSDArray<3>& bra,
               const btas::QSDArray<3>& ket)
{
  btas::QSDArray<4> escr;
  btas::QSDcontract(1.0, bra.conjugate(), shape(1), ket, shape(1), 1.0, escr);
  btas::QSDArray<4> eopr;
  btas::QSDpermute(escr, shape(0,2,1,3), eopr);
  btas::SDArray<1> sval;
  btas::QSDArray<3> lvec;
  btas::QSDArray<3> rvec;
  btas::QSDgesvd(btas::LeftArrow, eopr, sval, lvec, rvec, 0);
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\ttransfer operator :: ";
  cout << setprecision(8) << fixed << eopr << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
  cout << "\tsingular values   :: ";
  cout << setprecision(8) << fixed << sval << endl;
  cout << "\t----------------------------------------------------------------------------------------------------" << endl;
}



