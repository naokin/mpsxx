#include <iostream>
#include <iomanip>
using namespace std;

#include <legacy/Dblas.h>
#include <legacy/Dlapack.h>
#include <legacy/Dcontract.h>
#include <legacy/Ddiagonal.h>
#include "btas_template_specialize.h"

#include "driver.h"
using btas::shape;
using btas::NoTrans;
using btas::Trans;

void prototype::ComputeGuess
(bool forward, const btas::DArray<3>& mps0,
               const btas::DArray<3>& wfn0,
               const btas::DArray<3>& mps1,
                     btas::DArray<3>& wfn1)
{
  if(forward) {
    btas::DArray<2> lres;
    btas::Dgemm(  Trans, NoTrans, 1.0, mps0, wfn0, 1.0, lres);
    btas::Dgemm(NoTrans, NoTrans, 1.0, lres, mps1, 1.0, wfn1);
  }
  else {
    btas::DArray<2> rres;
    btas::Dgemm(NoTrans,   Trans, 1.0, wfn0, mps0, 1.0, rres);
    btas::Dgemm(NoTrans, NoTrans, 1.0, mps1, rres, 1.0, wfn1);
  }
}

void prototype::Canonicalize
(bool forward, const btas::DArray<3>& wfn0,
                     btas::DArray<3>& mps0, int M)
{
  btas::DArray<1> s;
  if(forward) {
    btas::DArray<2> v;
    Decompose(wfn0, s, mps0, v, M);
  }
  else {
    btas::DArray<2> u;
    Decompose(wfn0, s, u, mps0, M);
  }
  cout << "\tSingular values @ Canonicalization" << endl;
  for(btas::DArray<1>::iterator it = s.begin(); it != s.end(); ++it)
    cout << "\t\t\t" << setw(16) << setprecision(12) << fixed << *it << endl;
}

void prototype::Canonicalize
(bool forward, const btas::DArray<3>& wfn0,
                     btas::DArray<3>& mps0,
                     btas::DArray<3>& mps1)
{
  int n0 = wfn0.extent(0);
  int n1 = wfn0.extent(1);
  int n2 = wfn0.extent(2);

  btas::DArray<1> s;
  if(forward) {
    int ns = n2;
    int nt = n0*n1-ns;
    btas::DArray<3> u;
    btas::DArray<2> v;
    btas::Dgesvd(wfn0, s, u, v, ClapackCalcVector);
    mps0.resize(n0, n1, ns);
    for(int i = 0; i < n0; ++i)
      for(int j = 0; j < n1; ++j)
        for(int k = 0; k < ns; ++k)
          mps0(i, j, k) = u(i, j, k);
    mps1.resize(n0, n1, nt);
    for(int i = 0; i < n0; ++i)
      for(int j = 0; j < n1; ++j)
        for(int k = 0; k < nt; ++k)
          mps1(i, j, k) = u(i, j, ns+k);
  }
  else {
    int ns = n0;
    int nt = n1*n2-ns;
    btas::DArray<2> u;
    btas::DArray<3> v;
    btas::Dgesvd(wfn0, s, u, v, ClapackCalcVector);
    mps0.resize(ns, n1, n2);
    for(int i = 0; i < ns; ++i)
      for(int j = 0; j < n1; ++j)
        for(int k = 0; k < n2; ++k)
          mps0(i, j, k) = v(i, j, k);
    mps1.resize(nt, n1, n2);
    for(int i = 0; i < nt; ++i)
      for(int j = 0; j < n1; ++j)
        for(int k = 0; k < n2; ++k)
          mps1(i, j, k) = v(ns+i, j, k);
  }
  cout << "\tSingular values @ Canonicalization" << endl;
  for(btas::DArray<1>::iterator it = s.begin(); it != s.end(); ++it)
    cout << "\t\t\t" << setw(16) << setprecision(12) << fixed << *it << endl;
}

void prototype::Canonicalize
(bool forward, const btas::DArray<4>& wfnx,
                     btas::DArray<3>& mps0,
                     btas::DArray<3>& wfn1, int M)
{
  btas::DArray <1> s;
  if(forward) {
    Decompose(wfnx, s, mps0, wfn1, M);
    btas::Ddidm(s, wfn1);
  }
  else {
    Decompose(wfnx, s, wfn1, mps0, M);
    btas::Ddimd(wfn1, s);
  }
  cout << "\tSingular values @ Canonicalization" << endl;
  for(btas::DArray<1>::iterator it = s.begin(); it != s.end(); ++it)
    cout << "\t\t\t" << setw(16) << setprecision(12) << fixed << *it << endl;
}

void prototype::Renormalize
(bool forward, const btas::DArray<4>& mpo0,
               const btas::DArray<3>& opr0,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket0,
                     btas::DArray<3>& opr1)
{
  if(forward) {
    btas::DArray<4> scr1;
    btas::Dcontract(1.0, opr0, shape(0), bra0, shape(0), 1.0, scr1);
    btas::DArray<4> scr2;
    btas::Dcontract(1.0, scr1, shape(0, 2), mpo0, shape(0, 1), 1.0, scr2);
    btas::Dcontract(1.0, scr2, shape(0, 2), ket0, shape(0, 1), 1.0, opr1);
  }
  else {
    btas::DArray<4> scr1;
    btas::Dcontract(1.0, bra0, shape(2), opr0, shape(0), 1.0, scr1);
    btas::DArray<4> scr2;
    btas::Dcontract(1.0, scr1, shape(1, 2), mpo0, shape(1, 3), 1.0, scr2);
    btas::Dcontract(1.0, scr2, shape(3, 1), ket0, shape(1, 2), 1.0, opr1);
  }
}

void prototype::Renormalize
(bool forward, const btas::DArray<2>& opr0,
               const btas::DArray<3>& bra0,
               const btas::DArray<3>& ket0,
                     btas::DArray<2>& opr1)
{
  if(forward) {
    btas::DArray<3> scr1;
    btas::Dcontract(1.0, opr0, shape(0), bra0, shape(0), 1.0, scr1);
    btas::Dcontract(1.0, scr1, shape(0, 1), ket0, shape(0, 1), 1.0, opr1);
  }
  else {
    btas::DArray<3> scr1;
    btas::Dcontract(1.0, bra0, shape(2), opr0, shape(0), 1.0, scr1);
    btas::Dcontract(1.0, scr1, shape(1, 2), ket0, shape(1, 2), 1.0, opr1);
  }
}

void prototype::ComputeInverseGauge
(              const btas::DArray<2>& gauge_0,
                     btas::DArray<2>& gauge_i)
{
  btas::DArray<1> s;
  btas::DArray<2> u;
  btas::DArray<2> v;
  btas::Dgesvd(gauge_0, s, u, v);
  for(btas::DArray<1>::iterator it = s.begin(); it != s.end(); ++it) {
    if(fabs(*it) >= 1.0e-20)
      *it = 1.0/(*it);
    else
      *it = 0.0;
  }
  btas::Ddimd(u, s);
  gauge_i.free();
  btas::Dgemm(btas::Trans, btas::Trans, 1.0, v, u, 1.0, gauge_i);
}

void prototype::Normalize(btas::DArray<3>& wfn0)
{
  double norm = btas::Ddot(wfn0, wfn0);
  if(norm >= 1.0e-20) {
    btas::Dscal(1.0/sqrt(norm), wfn0);
  }
}

void prototype::Orthogonalize(const btas::DArray<3>& wfn0, btas::DArray<3>& wfn1)
{
  double overlap = btas::Ddot(wfn0, wfn1);
  Daxpy(-overlap, wfn0, wfn1);
}

void prototype::Orthogonalize(bool forward, const btas::DArray<3>& mps0, btas::DArray<3>& wfn1)
{
  if(forward) {
    btas::DArray<2> proj;
    btas::Dgemm(  Trans, NoTrans, 1.0, mps0, wfn1, 1.0, proj);
    btas::Dgemm(NoTrans, NoTrans,-1.0, mps0, proj, 1.0, wfn1);
  }
  else {
    btas::DArray<2> proj;
    btas::Dgemm(NoTrans,   Trans, 1.0, wfn1, mps0, 1.0, proj);
    btas::Dgemm(NoTrans, NoTrans,-1.0, proj, mps0, 1.0, wfn1);
  }
}

void prototype::ComputeDiagonal
(              const btas::DArray<4>& mpo0,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
                     btas::DArray<3>& diag)
{
  btas::DArray<3> mpo0_diag;
  btas::DArray<2> lopr_diag;
  btas::DArray<2> ropr_diag;

  btas::Ddiagonal(mpo0, shape(1, 2), mpo0_diag);
  btas::Ddiagonal(lopr, shape(0, 2), lopr_diag);
  btas::Ddiagonal(ropr, shape(0, 2), ropr_diag);

  btas::DArray<3> scr1;
  btas::Dcontract(1.0, lopr_diag, shape(1), mpo0_diag, shape(0), 1.0, scr1);
  btas::DArray<3> scr2;
  btas::Dcontract(1.0, scr1,      shape(2), ropr_diag, shape(1), 1.0, scr2);
  btas::Dcopy(scr2, diag);
}

void prototype::ComputeDiagonal
(              const btas::DArray<4>& lmpo,
               const btas::DArray<4>& rmpo,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
                     btas::DArray<4>& diag)
{
  btas::DArray<3> lmpo_diag;
  btas::DArray<3> rmpo_diag;
  btas::DArray<2> lopr_diag;
  btas::DArray<2> ropr_diag;

  btas::Ddiagonal(lmpo, shape(1, 2), lmpo_diag);
  btas::Ddiagonal(rmpo, shape(1, 2), rmpo_diag);
  btas::Ddiagonal(lopr, shape(0, 2), lopr_diag);
  btas::Ddiagonal(ropr, shape(0, 2), ropr_diag);

  btas::DArray<3> scr1;
  btas::Dcontract(1.0, lopr_diag, shape(1), lmpo_diag, shape(0), 1.0, scr1);
  btas::DArray<4> scr2;
  btas::Dcontract(1.0, scr1,      shape(2), rmpo_diag, shape(0), 1.0, scr2);
  btas::DArray<4> scr3;
  btas::Dcontract(1.0, scr2,      shape(3), ropr_diag, shape(1), 1.0, scr3);
  btas::Dcopy(scr3, diag);
}

void prototype::ComputeSigmaVector // compute: H | psi >
(              const btas::DArray<4>& mpo0,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
               const btas::DArray<3>& wfn0,
                     btas::DArray<3>& sgv0)
{
  btas::DArray<4> scr1;
  btas::Dcontract(1.0, lopr, shape(2), wfn0, shape(0), 1.0, scr1);
  btas::DArray<4> scr2;
  btas::Dcontract(1.0, scr1, shape(1, 2), mpo0, shape(0, 2), 1.0, scr2);
  btas::Dcontract(1.0, scr2, shape(3, 1), ropr, shape(1, 2), 1.0, sgv0);
}

void prototype::ComputeSigmaVectorConj // compute: < psi | H  ... but returned as ket state | sgv >
(              const btas::DArray<4>& mpo0,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
               const btas::DArray<3>& wfn0,
                     btas::DArray<3>& sgv0)
{
  btas::DArray<4> scr1;
  btas::Dcontract(1.0, lopr, shape(0), wfn0, shape(0), 1.0, scr1);
  btas::DArray<4> scr2;
  btas::Dcontract(1.0, scr1, shape(0, 2), mpo0, shape(0, 1), 1.0, scr2);
  btas::Dcontract(1.0, scr2, shape(1, 3), ropr, shape(0, 1), 1.0, sgv0);
}

void prototype::ComputeSigmaVector
(              const btas::DArray<4>& lmpo,
               const btas::DArray<4>& rmpo,
               const btas::DArray<3>& lopr,
               const btas::DArray<3>& ropr,
               const btas::DArray<4>& wfn0,
                     btas::DArray<4>& sgv0)
{
  btas::DArray<5> scr1;
  btas::Dcontract(1.0, lopr, shape(2), wfn0, shape(0), 1.0, scr1);
  btas::DArray<5> scr2;
  btas::Dcontract(1.0, scr1, shape(1, 2), lmpo, shape(0, 2), 1.0, scr2);
  btas::DArray<5> scr3;
  btas::Dcontract(1.0, scr2, shape(4, 1), rmpo, shape(0, 2), 1.0, scr3);
  btas::Dcontract(1.0, scr3, shape(4, 1), ropr, shape(1, 2), 1.0, sgv0);
}

int prototype::ComputeOrthogonalTransform
(              const btas::DArray<2>& s_matrix,
                     btas::DArray<2>& u_matrix)
{
  int nrows = s_matrix.rows();
  int ncols = s_matrix.cols();
  // diagonalize s_subspace
  btas::DArray<1> d;
  btas::DArray<2> p_matrix;
  btas::Dsyev(s_matrix, d, p_matrix);

  cout << "\tprinting eigenvalues of overlap matrix" << endl;
  for(int i = 0; i < d.size(); ++i)
    cout << "\t\tState [ " << setw(2) << i + 1 << " ]: " << setw(16) << fixed << setprecision(12) << d(i) << endl;

  // check positive eigenvalues;
//int nimag = 0;
//for(; nimag < d.size(); ++nimag) if(d(nimag) >= 1.0e-20) break;
  int nmin = 0;
  for(; nmin < d.size(); ++nmin)
    if(d(nmin) >= 1.0e-16) break;
  int nmax = ncols-1-nmin;
//int nmax = ncols-1;
//for(; nmax >= 0; --nmax)
//  if(d(nmax) < 1.0e+4) break;
  ncols = nmax - nmin + 1;

  u_matrix.resize(nrows, ncols);

  int icol = 0;
  for(int i = nmax; i >= nmin; --i, ++icol) {
    double di = 1.0/sqrt(d(i));
    for(int j = 0; j < nrows; ++j)
      u_matrix(j, icol) = p_matrix(i, j)*di;
  }

//int icol = 0;
//for(int i = nimag; i < p_matrix.rows(); ++i, ++icol) {
//  di(icol) = 1.0/sqrt(d(i));
//  for(int j = 0; j < p_matrix.cols(); ++j) {
//    u_matrix(j, icol) = p_matrix(i, j);
//  }
//}
//btas::Ddimd(u_matrix, di);

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
  alpha.free();
  btas::Dgemm(NoTrans, Trans, 1.0, u_matrix, t_matrix, 1.0, alpha);

  return ncols;
}

int prototype::ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm)
{
  int nrows = a_subspace.rows();
  int ncols = a_subspace.cols();

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
//cout << "debug: u_matrix: [" << u_matrix.rows() << "," << u_matrix.cols() << "]" << endl;
//for(int i = 0; i < u_matrix.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < u_matrix.cols(); ++j) {
//    cout << setw(8) << fixed << u_matrix(i, j);
//  }
//  cout << endl;
//}

  btas::DArray<2> s_reduce;
  TransformedMatrix(s_matrix, u_matrix, s_reduce);

  btas::DArray<1> d;
  btas::DArray<2> p;
  btas::Dsyev(s_reduce, d, p); // h_matrix as metric

  // check positive eigenvalues;
  int nmin = 0;
  for(; nmin < d.size(); ++nmin)
    if(d(nmin) >= 1.0e-20) break;
  int nmax = ncols-1;
  for(; nmax >= 0; --nmax)
    if(d(nmax) < 1.0e+20) break;
  ncols = nmax - nmin + 1;

  cout << "\tprinting eigenvalues (1/h)" << endl;
  for(int i = 0; i < d.size(); ++i)
  cout << "\t\t#" << setw(2) << i << ": " << setprecision(8) << setw(16) << scientific << d(i) << endl;
//cout << "debug: (after diagonalize) nmin = " << nmin << ", nmax = " << nmax << " ncols = " << ncols << endl;
  btas::DArray<2> x_matrix;
  btas::Dgemm(NoTrans, Trans, 1.0, u_matrix, p, 1.0, x_matrix);

// DEBUG
//cout.precision(4);
//cout << "debug: x_matrix: [" << x_matrix.rows() << "," << x_matrix.cols() << "]" << endl;
//for(int i = 0; i < x_matrix.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < x_matrix.cols(); ++j) {
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
//cout << "debug: alpha: [" << alpha.rows() << "," << alpha.cols() << "]" << endl;
//for(int i = 0; i < alpha.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < alpha.cols(); ++j) {
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
//for(int i = 0; i < h_diag.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < h_diag.cols(); ++j) {
//    cout << setw(8) << fixed << h_diag(i, j);
//  }
//  cout << endl;
//}
//btas::DArray<2> s_diag;
//btas::DArray<2> s_diagtmp;
//btas::Dgemm(NoTrans, NoTrans, 1.0, s_matrix, alpha, 1.0, s_diagtmp);
//btas::Dgemm(  Trans, NoTrans, 1.0, alpha, s_diagtmp, 1.0, s_diag);
//cout << "DEBUG: S(diagonalized)" << endl;
//for(int i = 0; i < s_diag.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < s_diag.cols(); ++j) {
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
//cout << "debug: alphaRe: [" << alphaRe.rows() << "," << alphaRe.cols() << "]" << endl;
//for(int i = 0; i < alphaRe.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < alphaRe.cols(); ++j) {
//    cout << setw(8) << fixed << alphaRe(i, j);
//  }
//  cout << endl;
//}
//cout << "debug: alphaIm: [" << alphaIm.rows() << "," << alphaIm.cols() << "]" << endl;
//for(int i = 0; i < alphaIm.rows(); ++i) {
//  cout << "\t";
//  for(int j = 0; j < alphaIm.cols(); ++j) {
//    cout << setw(8) << fixed << alphaIm(i, j);
//  }
//  cout << endl;
//}

  return ncols;
}




//
// FOR TEST
//

int prototype::TEST::ComputeOrthogonalTransform
(              const btas::DArray<2>& s_matrix,
                     btas::DArray<2>& u_matrix)
{
  int nrows = s_matrix.rows();
  int ncols = s_matrix.cols();
  // diagonalize s_subspace
  btas::DArray<1> d;
  btas::DArray<2> p_matrix;
  btas::Dsyev(s_matrix, d, p_matrix);

  cout << "\tprinting eigenvalues of overlap matrix" << endl;
  std::vector<int> nnz_index;
  for(int i = 0; i < d.size(); ++i) {
    if(fabs(d(i)) > 1.0e-16) nnz_index.push_back(i);
    cout << "\t\tState [ " << setw(2) << i + 1 << " ]: " << setw(16) << fixed << setprecision(12) << d(i) << endl;
  }

  ncols = nnz_index.size();

  u_matrix.resize(nrows, ncols);
  for(int icol = ncols-1; icol >= 0; --icol) {
    int i = nnz_index[icol];
    double di = 1.0/sqrt(fabs(d(i)));
    for(int j = 0; j < nrows; ++j)
      u_matrix(j, icol) = p_matrix(i, j)*di;
  }

  return ncols;
}

int prototype::TEST::ComputeEigenvalues
(              const btas::DArray<2>& a_subspace,
               const btas::DArray<2>& b_subspace,
               const btas::DArray<2>& s_subspace,
               const btas::DArray<2>& d_subspace,
                     btas::DArray<1>& eigvs,
                     btas::DArray<2>& alphaRe,
                     btas::DArray<2>& alphaIm)
{
  int nrows = a_subspace.rows();
  int ncols = a_subspace.cols();

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
  ncols = prototype::TEST::ComputeOrthogonalTransform(s_matrix, u_matrix);

  btas::DArray<2> h_reduce;
  TransformedMatrix(h_matrix, u_matrix, h_reduce);

  btas::DArray<2> s_reduce;
  TransformedMatrix(s_matrix, u_matrix, s_reduce);

  btas::DArray<1> d;
  btas::DArray<2> p;
  btas::Dsyev(h_reduce, d, p);

  btas::DArray<2> s_tmp;
  btas::Dgemm(NoTrans, NoTrans, 1.0, p, s_reduce, 1.0, s_tmp);
  btas::DArray<2> s;
  btas::Dgemm(NoTrans,   Trans, 1.0, s_tmp, p, 1.0, s);

  cout << "\t***** test computation *****" << endl;
  cout << "\tprinting eigenvalues" << endl;
  for(int i = 0; i < d.size(); ++i) {
    cout << "\t" << setw(12) << setprecision(8) << fixed << d(i) << endl;
  }
  cout << "\tprinting reduced metric" << endl;
  for(int i = 0; i < s.rows(); ++i) {
    cout << "\t";
    for(int j = 0; j < s.cols(); ++j) {
      cout << setw(12) << setprecision(8) << fixed << s(i, j);
    }
    cout << endl;
  }

//// check positive eigenvalues;
//int nmin = 0;
//for(; nmin < d.size(); ++nmin)
//  if(d(nmin) >= 1.0e-16) break;
//ncols -= nmin;

//cout << "\tprinting eigenvalues (1/h)" << endl;
//for(int i = 0; i < d.size(); ++i)
//cout << "\t\t#" << setw(2) << i << ": " << setprecision(8) << setw(16) << scientific << d(i) << endl;
//btas::DArray<2> x_matrix;
//btas::Dgemm(NoTrans, Trans, 1.0, u_matrix, p, 1.0, x_matrix);

//ncols /= 2;
//btas::DArray<2> alpha(2*nrows, ncols);
//eigvs.resize(ncols);

//int icol = 0;
//for(int i = nmin/2; i < ncols; ++i, ++icol) {
//  eigvs(icol) = d(2*i);
//  for(int j = 0; j < 2*nrows; ++j)
//    alpha(j, icol) = x_matrix(j, 2*i);
//}

//alphaRe.resize(nrows, ncols);
//alphaIm.resize(nrows, ncols);
//for(int i = 0; i < nrows; ++i) {
//  for(int j = 0; j < ncols; ++j) {
//    alphaRe(i, j) = alpha(i,       j);
//    alphaIm(i, j) = alpha(i+nrows, j);
//  }
//}

  return ncols;
}



