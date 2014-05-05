#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <dense.h>
using namespace std;
using namespace btas;

#define DEBUG(msg) { cout << "DEBUG:" << msg << endl; }

template<int N>
void normalize(DTensor<N>& t)
{
  double norm = Ddot(t, t);
  Dscal(1.0/sqrt(norm),t);
}

void decompose(const bool& forward, const DTensor<4>& psi, DTensor<3>& lmps, DTensor<3>& rmps, const int& M)
{
  DTensor<1> stmp;
  DTensor<3> ltmp;
  DTensor<3> rtmp;
  Dgesvd(psi,shape(2,3),stmp,ltmp,rtmp);
  int M0 = stmp.size();
  for(int m = 0; m < stmp.size(); ++m)
    if(stmp(m) < 1.0e-12) {
      M0 = m; break;
    }
  if(M > 0) M0 = min(M0,M);
  int Ml = ltmp.extent(0);
  int dl = ltmp.extent(1);
  int dr = rtmp.extent(1);
  int Mr = rtmp.extent(2);
  // truncated singular values
  DTensor<1> sval(M0);
  for(int i = 0; i < M0; ++i) sval(i) = stmp(i);
  normalize(sval);
  // truncated left mps
  lmps.resize(Ml,dl,M0);
  for(int i = 0; i < Ml; ++i)
    for(int j = 0; j < dl; ++j)
      for(int k = 0; k < M0; ++k)
        lmps(i,j,k) = ltmp(i,j,k);
  // truncated right mps
  rmps.resize(M0,dr,Mr);
  for(int i = 0; i < M0; ++i)
    for(int j = 0; j < dr; ++j)
      for(int k = 0; k < Mr; ++k)
        rmps(i,j,k) = rtmp(i,j,k);
  if(forward)
    Dright_update(sval,rmps);
  else
    Dleft_update(lmps,sval);
}

void renormalize(const bool& forward, const DTensor<4>& mpo, const DTensor<3>& mps, DTensor<3>& storage)
{
  if(forward) {
    DTensor<4> scr1;
    Dcontract(1.0,storage,shape(0),mps,shape(0),1.0,scr1);
    DTensor<4> scr2;
    Dcontract(1.0,scr1,shape(0,2),mpo,shape(0,2),1.0,scr2);
    storage.free();
    Dcontract(1.0,scr2,shape(0,3),mps,shape(0,1),1.0,storage);
  }
  else {
    DTensor<4> scr1;
    Dcontract(1.0,mps,shape(2),storage,shape(0),1.0,scr1);
    DTensor<4> scr2;
    Dcontract(1.0,scr1,shape(1,2),mpo,shape(2,1),1.0,scr2);
    storage.free();
    Dcontract(1.0,scr2,shape(3,1),mps,shape(1,2),1.0,storage);
  }
}

void renormalize(const bool& forward, const DTensor<3>& mps, DTensor<2>& storage)
{
  if(forward) {
    DTensor<3> scr1;
    Dcontract(1.0,storage,shape(0),mps,shape(0),1.0,scr1);
    storage.free();
    Dcontract(1.0,scr1,shape(0,1),mps,shape(0,1),1.0,storage);
  }
  else {
    DTensor<3> scr1;
    Dcontract(1.0,mps,shape(2),storage,shape(0),1.0,scr1);
    storage.free();
    Dcontract(1.0,scr1,shape(1,2),mps,shape(1,2),1.0,storage);
  }
}

double compute_energy(const bool& forward, const vector< DTensor<3> >& mpsts, const vector< DTensor<4> >& mpops)
{
  // suppose mpsts is right-canonical form
  int len = mpsts.size();
  DTensor<3> sops(1,1,1); sops = 1.0;
  DTensor<2> snrm(1,1);   snrm = 1.0;
  if(forward)
    for(int i = 0; i < len; ++i) {
      renormalize(forward,mpops[i],mpsts[i],sops);
      renormalize(forward,mpsts[i],snrm);
    }
  else
    for(int i = len-1; i >= 0; --i) {
      renormalize(forward,mpops[i],mpsts[i],sops);
      renormalize(forward,mpsts[i],snrm);
    }
  DTensor<3> cops(1,1,1); cops = 1.0;
  DTensor<2> cnrm(1,1);   cnrm = 1.0;
  double en = Ddot(sops,cops);
  double nm = Ddot(snrm,cnrm);
  return en/nm;
}

void tdmrg_evolving(const bool& forward, const DTensor<4>& ops, DTensor<3>& lmps, DTensor<3>& rmps, const int& M)
{
  DTensor<4> psi;
  Dcontract(1.0,lmps,shape(2),rmps,shape(0),1.0,psi);
  DTensor<4> phi;
  {
    DTensor<4> scr;
    Dcontract(1.0,psi,shape(1,2),ops,shape(2,3),1.0,scr);
    Dpermute(scr,shape(0,2,3,1),phi);
  }
  decompose(forward,phi,lmps,rmps,M);
}

void tdmrg_skipping(const bool& forward, DTensor<3>& lmps, DTensor<3>& rmps, const int& M)
{
  DTensor<4> psi;
  Dcontract(1.0,lmps,shape(2),rmps,shape(0),1.0,psi);
  decompose(forward,psi,lmps,rmps,M);
}

void tdmrg_sweeping(bool& forward, const DTensor<4>& hodd, const DTensor<4>& heven, vector< DTensor<3> >& mpsts, const int& M)
{
  int len = mpsts.size();
  int n; bool finalize;
  if(len % 2 == 0) {
    n = len/2 - 1;
    finalize = true;
  }
  else {
    n = len/2;
    finalize = false;
  }
  if(forward) {
    // evolving odd # sites
    for(int i = 0; i < n; ++i) {
      int i2 = 2 * i;
      tdmrg_evolving(forward, hodd,  mpsts[i2],   mpsts[i2+1], M);
      tdmrg_skipping(forward,        mpsts[i2+1], mpsts[i2+2], M);
    }
    if(finalize) {
      int i2 = 2 * n;
      tdmrg_evolving(forward, hodd,  mpsts[i2],   mpsts[i2+1], M);
      forward = !forward;
      tdmrg_skipping(forward,        mpsts[i2],   mpsts[i2+1], M);
    }
    // evolving even # sites
    for(int i = n; i > 0; --i) {
      int i2 = 2 * i;
      tdmrg_evolving(forward, heven, mpsts[i2-1], mpsts[i2],   M);
      tdmrg_skipping(forward,        mpsts[i2-2], mpsts[i2-1], M);
    }
    forward = !forward;
    // evolving odd # sites
    for(int i = 0; i < n; ++i) {
      int i2 = 2 * i;
      tdmrg_evolving(forward, hodd,  mpsts[i2],   mpsts[i2+1], M);
      tdmrg_skipping(forward,        mpsts[i2+1], mpsts[i2+2], M);
    }
    if(finalize) {
      int i2 = 2 * n;
      tdmrg_evolving(forward, hodd,  mpsts[i2],   mpsts[i2+1], M);
    }
    forward = !forward;
  }
  else {
    // evolving odd # sites
    if(finalize) {
      int i2 = 2 * n;
      tdmrg_evolving(forward, hodd,  mpsts[i2],   mpsts[i2+1], M);
    }
    for(int i = n; i > 0; --i) {
      int i2 = 2 * i;
      tdmrg_skipping(forward,        mpsts[i2-1], mpsts[i2],   M);
      tdmrg_evolving(forward, hodd,  mpsts[i2-2], mpsts[i2-1], M);
    }
    forward = !forward;
    // evolving even # sites
    for(int i = 0; i < n; ++i) {
      int i2 = 2 * i;
      tdmrg_skipping(forward,        mpsts[i2],   mpsts[i2+1], M);
      tdmrg_evolving(forward, heven, mpsts[i2+1], mpsts[i2+2], M);
    }
    if(finalize) {
      int i2 = 2 * n;
      tdmrg_skipping(forward,        mpsts[i2],   mpsts[i2+1], M);
      forward = !forward;
      tdmrg_evolving(forward, hodd,  mpsts[i2],   mpsts[i2+1], M);
    }
    // evolving odd # sites
    for(int i = n; i > 0; --i) {
      int i2 = 2 * i;
      tdmrg_skipping(forward,        mpsts[i2-1], mpsts[i2],   M);
      tdmrg_evolving(forward, hodd,  mpsts[i2-2], mpsts[i2-1], M);
    }
    forward = !forward;
  }
}

int imag_tdmrg(ostream& fout, vector<DTensor<3> >& mpsts, bool restart,
               int L, int M, double J, double Jz, int N, double dt)
{
  fout.setf(ios::fixed,ios::floatfield);
  fout.precision(4);
  fout << "====================  LATTICE INFO  ====================" << endl;
  fout << "\tL =" << setw(4) << L << ", M ="  << setw(4) << M
       << ", J =" << setw(7) << J << ", Jz =" << setw(7) << Jz << endl;

  DTensor<4> hodd(2,2,2,2);
  hodd = 0.0;
  {
    double Xz = exp(-Jz * 0.5 * dt / 4);
    double Yz = exp( Jz * 0.5 * dt / 4);
    double X  = exp(-J  * 0.5 * dt / 2);
    double Y  = exp( J  * 0.5 * dt / 2);
    hodd(0,0,0,0) = Xz;
    hodd(0,1,0,1) = Yz * (X + Y) / 2;
    hodd(0,1,1,0) = Yz * (X - Y) / 2;
    hodd(1,0,0,1) = Yz * (X - Y) / 2;
    hodd(1,0,1,0) = Yz * (X + Y) / 2;
    hodd(1,1,1,1) = Xz;
//  hodd(0,0,0,0) = 1.0 - Jz * 0.5 * dt / 4;
//  hodd(0,1,0,1) = 1.0 + Jz * 0.5 * dt / 4;
//  hodd(0,1,1,0) =     - J  * 0.5 * dt / 2;
//  hodd(1,0,0,1) =     - J  * 0.5 * dt / 2;
//  hodd(1,0,1,0) = 1.0 + Jz * 0.5 * dt / 4;
//  hodd(1,1,1,1) = 1.0 - Jz * 0.5 * dt / 4;
  }

  DTensor<4> heven(2,2,2,2);
  heven = 0.0;
  {
    double Xz = exp(-Jz * dt / 4);
    double Yz = exp( Jz * dt / 4);
    double X  = exp(-J  * dt / 2);
    double Y  = exp( J  * dt / 2);
    heven(0,0,0,0) = Xz;
    heven(0,1,0,1) = Yz * (X + Y) / 2;
    heven(0,1,1,0) = Yz * (X - Y) / 2;
    heven(1,0,0,1) = Yz * (X - Y) / 2;
    heven(1,0,1,0) = Yz * (X + Y) / 2;
    heven(1,1,1,1) = Xz;
//  heven(0,0,0,0) = 1.0 - Jz * dt / 4;
//  heven(0,1,0,1) = 1.0 + Jz * dt / 4;
//  heven(0,1,1,0) =     - J  * dt / 2;
//  heven(1,0,0,1) =     - J  * dt / 2;
//  heven(1,0,1,0) = 1.0 + Jz * dt / 4;
//  heven(1,1,1,1) = 1.0 - Jz * dt / 4;
  }

  // initialize mpstates
  // create anti-parallel spins state : ababab...
  // and this will always be kept as right-canonical
  if(!restart) {
    int M0 = 1;
    mpsts.resize(L);
    mpsts[0].resize(1,2,M0);
    mpsts[0](0,0,0) = 1.0;
    for(int i = 1; i < L-1; ++i) {
      mpsts[i].resize(M0,2,M0);
      if(i & 1) mpsts[i](0,1,0) = 1.0;
      else      mpsts[i](0,0,0) = 1.0;
    }
    mpsts[L-1].resize(M0,2,1);
    if(L & 1) mpsts[L-1](0,0,0) = 1.0;
    else      mpsts[L-1](0,1,0) = 1.0;
  }
  // mpo's
  vector< DTensor<4> > mpops;
  mpops.resize(L);
  {
    // [ 0  JS+  JS-  JzSz  I ]
    DTensor<4> mpo0(1,5,2,2);
    mpo0 = 0.0;
    mpo0(0,1,0,1) = J /2; // S+
    mpo0(0,2,1,0) = J /2; // S-
    mpo0(0,3,0,0) = Jz/2; // Sz(a)
    mpo0(0,3,1,1) =-Jz/2; // Sz(b)
    mpo0(0,4,0,0) = 1.00; // I(a)
    mpo0(0,4,1,1) = 1.00; // I(b)
    Dcopy(mpo0,mpops[0]);
  }
  for(int i = 1; i < L-1; ++i) {
    // [ I   0    0    0    0 ]
    // [ S-  0    0    0    0 ]
    // [ S+  0    0    0    0 ]
    // [ Sz  0    0    0    0 ]
    // [ 0  JS+  JS-  JzSz  I ]
    DTensor<4> mpoI(5,5,2,2);
    mpoI = 0.0;
    mpoI(0,0,0,0) = 1.00; // I(a)
    mpoI(0,0,1,1) = 1.00; // I(b)
    mpoI(1,0,1,0) = 1.00; // S-
    mpoI(2,0,0,1) = 1.00; // S+
    mpoI(3,0,0,0) = 0.50; // Sz(a)
    mpoI(3,0,1,1) =-0.50; // Sz(b)
    mpoI(4,1,0,1) = J /2; // S+
    mpoI(4,2,1,0) = J /2; // S-
    mpoI(4,3,0,0) = Jz/2; // Sz(a)
    mpoI(4,3,1,1) =-Jz/2; // Sz(b)
    mpoI(4,4,0,0) = 1.00; // I(a)
    mpoI(4,4,1,1) = 1.00; // I(b)
    Dcopy(mpoI,mpops[i]);
  }
  {
    // [ I  ]
    // [ S- ]
    // [ S+ ]
    // [ Sz ]
    // [ 0  ]
    DTensor<4> mpoN(5,1,2,2);
    mpoN = 0.0;
    mpoN(0,0,0,0) = 1.00; // I(a)
    mpoN(0,0,1,1) = 1.00; // I(b)
    mpoN(1,0,1,0) = 1.00; // S-
    mpoN(2,0,0,1) = 1.00; // S+
    mpoN(3,0,0,0) = 0.50; // Sz(a)
    mpoN(3,0,1,1) =-0.50; // Sz(b)
    Dcopy(mpoN,mpops[L-1]);
  }

  // imaginary time-evolution
  fout << "==================== TIME-EVOLUTION ====================" << endl;
  double e   = 0.0;
  double t   = 0.0;
  int    itr = 0;
  bool   fwd = true;
  while(itr < N) {
    if(itr % 100 == 0) {
      e = compute_energy(fwd,mpsts,mpops);
      fout.precision(4);
      fout << "\tt = " << setw(8) << t;
      fout.precision(16);
      fout << " : E = " << setw(24) << e << endl;
    }
    // time-evolution
    tdmrg_sweeping(fwd,hodd,heven,mpsts,M);
    t += dt; ++itr;
  }
  e = compute_energy(fwd,mpsts,mpops);
  fout.precision(4);
  fout << "\tt = " << setw(8) << t;
  fout.precision(16);
  fout << " : E = " << setw(24) << e << endl;
  fout << "==================== FINISHED STEPS ====================" << endl;

  return 0;
}
