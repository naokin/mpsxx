#include <vector>
#include <cmath>
#include <legacy/Dblas.h>
#include <legacy/Dlapack.h>
#include <legacy/Dcontract.h>
#include "driver.h"
using namespace std;
using namespace btas;

// canonicalize / 1-site algorithm
double canonicalize(bool forward, const btas::DArray<3>& wfn0,
                                        btas::DArray<3>& mps0,
                                        btas::DArray<3>& wfn1, int M)
{
  if(forward) {
    DArray<1> s;
    DArray<2> v;
    Dgesvd(wfn0, s, mps0, v);
    Dright_update(s, v);
    DArray<3> wfnx;
    Dgemm(NoTrans, NoTrans, 1.0, v, wfn1, 1.0, wfnx);
    Dcopy(wfnx, wfn1);
  }
  else {
    DArray<1> s;
    DArray<2> u;
    Dgesvd(wfn0, s, u, mps0);
    Dleft_update(u, s);
    DArray<3> wfnx;
    Dgemm(NoTrans, NoTrans, 1.0, wfn1, u, 1.0, wfnx);
    Dcopy(wfnx, wfn1);
  }
  return 0.0;
}

// canonicalize / 1-site algorithm / from system density matrix
double canonicalize(bool forward, const btas::DArray<4>& sdm0,
                                        btas::DArray<1>& val0,
                                        btas::DArray<3>& mps0,
                                        btas::DArray<3>& nul0, int M)
{
  DArray<2> valx;
  DArray<4> rotx;
  Dsyev(sdm0, valx, rotx);

  int ni = sdm0.extent(0);
  int nj = sdm0.extent(1);
  int M0 = valx.size();

  val0.resize(M0);
  DArray<2>::iterator ivx = valx.begin();
  for(int i = M0-1; i >= 0; --i, ++ivx) val0(i) = *ivx;

  int Mx = M0;
  for(int i = 0; i < M0; ++i)
    if(fabs(val0(i)) < 1.0e-12) {
      Mx = i;
      break;
    }

  if(M > 0 && Mx > M) Mx = M;

  int Mn = M0 - Mx; if(Mn < 0) Mn = 0;

  if(forward) {
    DArray<4>::iterator irx = rotx.begin();
    nul0.resize(ni, nj, Mn);
    for(int k = Mn-1; k >= 0; --k)
      for(int i = 0; i < ni; ++i)
        for(int j = 0; j < nj; ++j)
          nul0(i, j, k) = *(irx++);
    mps0.resize(ni, nj, Mx);
    for(int k = Mx-1; k >= 0; --k)
      for(int i = 0; i < ni; ++i)
        for(int j = 0; j < nj; ++j)
          mps0(i, j, k) = *(irx++);
  }
  else {
    DArray<4>::iterator irx = rotx.begin();
    nul0.resize(Mn, ni, nj);
    for(int k = Mn-1; k >= 0; --k)
      for(int i = 0; i < ni; ++i)
        for(int j = 0; j < nj; ++j)
          nul0(k, i, j) = *(irx++);
    mps0.resize(Mx, ni, nj);
    for(int k = Mx-1; k >= 0; --k)
      for(int i = 0; i < ni; ++i)
        for(int j = 0; j < nj; ++j)
          mps0(k, i, j) = *(irx++);
  }

  double weight = 0.0;
  for(int i = Mx; i < M0; ++i) weight += val0(i);

  return fabs(weight);
}

// renormalize
void   renormalize (bool forward, const btas::DArray<4>& mpo0,
                                  const btas::DArray<3>& str0,
                                  const btas::DArray<3>& bra0,
                                  const btas::DArray<3>& ket0,
                                        btas::DArray<3>& str1)
{
  if(forward) {
    DArray<4> scr1;
    Dcontract(1.0, str0, shape(0),    bra0, shape(0),    1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(0, 2), mpo0, shape(0, 1), 1.0, scr2);
    Dcontract(1.0, scr2, shape(0, 2), ket0, shape(0, 1), 1.0, str1);
  }
  else {
    DArray<4> scr1;
    Dcontract(1.0, bra0, shape(2),    str0, shape(0),    1.0, scr1);
    DArray<4> scr2;
    Dcontract(1.0, scr1, shape(1, 2), mpo0, shape(1, 3), 1.0, scr2);
    Dcontract(1.0, scr2, shape(3, 1), ket0, shape(1, 2), 1.0, str1);
  }
}

void   renormalize (bool forward, const btas::DArray<2>& str0,
                                  const btas::DArray<3>& bra0,
                                  const btas::DArray<3>& ket0,
                                        btas::DArray<2>& str1)
{
  if(forward) {
    DArray<3> scr1;
    Dcontract(1.0, str0, shape(0),    bra0, shape(0),    1.0, scr1);
    Dcontract(1.0, scr1, shape(0, 1), ket0, shape(0, 1), 1.0, str1);
  }
  else {
    DArray<3> scr1;
    Dcontract(1.0, bra0, shape(2),    str0, shape(0),    1.0, scr1);
    Dcontract(1.0, scr1, shape(1, 2), ket0, shape(1, 2), 1.0, str1);
  }
}

// compute diagonal elements
void   compute_h_diagonal        (const btas::DArray<4>& mpo0,
                                  const btas::DArray<3>& lstr,
                                  const btas::DArray<3>& rstr,
                                        btas::DArray<3>& diag)
{
  int Dl = mpo0.extent(0);
  int Dr = mpo0.extent(1);
  int d0 = mpo0.extent(2);
  int Ml = lstr.extent(0);
  int Mr = rstr.extent(0);
  DArray<3> mpo0_diag(Dl, d0, Dr);
  for(int i = 0; i < Dl; ++i)
    for(int j = 0; j < d0; ++j)
      for(int k = 0; k < Dr; ++k)
        mpo0_diag(i, j, k) = mpo0(i, j, j, k);
  DArray<2> lstr_diag(Ml, Dl);
  for(int i = 0; i < Ml; ++i)
    for(int j = 0; j < Dl; ++j)
      lstr_diag(i, j) = lstr(i, j, i);
  DArray<2> rstr_diag(Mr,Dr);
  for(int i = 0; i < Mr; ++i)
    for(int j = 0; j < Dr; ++j)
      rstr_diag(i, j) = rstr(i, j, i);
  DArray<3> scr;
  Dcontract(1.0, lstr_diag, shape(1), mpo0_diag, shape(0), 1.0, scr);
  diag.free();
  Dcontract(1.0, scr, shape(2), rstr_diag, shape(1), 1.0, diag);
}

// compute sigma vector
void   compute_sigma_vector      (const btas::DArray<4>& mpo0,
                                  const btas::DArray<3>& lstr,
                                  const btas::DArray<3>& rstr,
                                  const btas::DArray<3>& mps0,
                                        btas::DArray<3>& sgv0)
{
  DArray<4> scr1;
  Dcontract(1.0, lstr, shape(2),    mps0, shape(0),    1.0, scr1);
  DArray<4> scr2;
  Dcontract(1.0, scr1, shape(1, 2), mpo0, shape(0, 2), 1.0, scr2);
  Dcontract(1.0, scr2, shape(3, 1), rstr, shape(1, 2), 1.0, sgv0);
}

