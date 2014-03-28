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

template<int N>
void normalize(ZTensor<N>& t)
{
  Znormalize(t);
}

template<int N>
void orthogonalize(const ZTensor<N>& base, ZTensor<N>& sub)
{
  Zorthogonalize(base, sub);
}

void decompose(const bool& forward, const ZTensor<4>& psi, ZTensor<3>& lmps, ZTensor<3>& rmps, const int& M)
{
  DTensor<1> stmp;
  ZTensor<3> ltmp;
  ZTensor<3> rtmp;
  Zgesvd(psi,shape(2,3),stmp,ltmp,rtmp);
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
    Zright_update(sval,rmps);
  else
    Zleft_update(lmps,sval);
}

void renormalize(const bool& forward, const ZTensor<4>& mpo, const ZTensor<3>& mps, ZTensor<3>& storage)
{
  if(forward) {
    ZTensor<4> scr1;
    Zcontract(1.0,storage,shape(0),mps.conjugate(),shape(0),1.0,scr1);
    ZTensor<4> scr2;
    Zcontract(1.0,scr1,shape(0,2),mpo,shape(0,2),1.0,scr2);
    storage.free();
    Zcontract(1.0,scr2,shape(0,3),mps,shape(0,1),1.0,storage);
  }
  else {
    ZTensor<4> scr1;
    Zcontract(1.0,mps.conjugate(),shape(2),storage,shape(0),1.0,scr1);
    ZTensor<4> scr2;
    Zcontract(1.0,scr1,shape(1,2),mpo,shape(2,1),1.0,scr2);
    storage.free();
    Zcontract(1.0,scr2,shape(3,1),mps,shape(1,2),1.0,storage);
  }
}

void renormalize(const bool& forward, const ZTensor<3>& mps, ZTensor<2>& storage)
{
  if(forward) {
    ZTensor<3> scr1;
    Zcontract(1.0,storage,shape(0),mps.conjugate(),shape(0),1.0,scr1);
    storage.free();
    Zcontract(1.0,scr1,shape(0,1),mps,shape(0,1),1.0,storage);
  }
  else {
    ZTensor<3> scr1;
    Zcontract(1.0,mps.conjugate(),shape(2),storage,shape(0),1.0,scr1);
    storage.free();
    Zcontract(1.0,scr1,shape(1,2),mps,shape(1,2),1.0,storage);
  }
}

double compute_energy(const bool& forward, const vector< ZTensor<3> >& mpsts, const vector< ZTensor<4> >& mpops)
{
  // suppose mpsts is right-canonical form
  int len = mpsts.size();
  ZTensor<3> sops(1,1,1); sops = Complx(1.0, 0.0);
  ZTensor<2> snrm(1,1);   snrm = Complx(1.0, 0.0);
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
  ZTensor<3> cops(1,1,1); cops = Complx(1.0, 0.0);
  ZTensor<2> cnrm(1,1);   cnrm = Complx(1.0, 0.0);
  Complx en = Zdotc(sops,cops);
  Complx nm = Zdotc(snrm,cnrm);
  return en.real()/nm.real();
}

void tdmrg_evolving(const bool& forward, const ZTensor<4>& ops, ZTensor<3>& lmps, ZTensor<3>& rmps, const int& M)
{
  ZTensor<4> psi;
  Zcontract(1.0,lmps,shape(2),rmps,shape(0),1.0,psi);
  ZTensor<4> phi;
  {
    ZTensor<4> scr;
    Zcontract(1.0,psi,shape(1,2),ops,shape(2,3),1.0,scr);
    Zpermute(scr,shape(0,2,3,1),phi);
  }
  decompose(forward,phi,lmps,rmps,M);
}

void tdmrg_skipping(const bool& forward, ZTensor<3>& lmps, ZTensor<3>& rmps, const int& M)
{
  ZTensor<4> psi;
  Zcontract(1.0,lmps,shape(2),rmps,shape(0),1.0,psi);
  decompose(forward,psi,lmps,rmps,M);
}

void tdmrg_sweeping(bool& forward, const ZTensor<4>& hodd, const ZTensor<4>& heven, vector< ZTensor<3> >& mpsts, const int& M)
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
    // what is n ?
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

ZTensor<4> propagator(double J, double Jz, double h, double f, double t, double dt)
{
  //
  // compute propagator exp(-iHdt) :
  //
  // H = sum{<i, j>} [ J ( Si+Sj- + Si-Sj+ ) / 2 + Jz SizSjz ] + sum{i} h cos(f * t) Siz
  //
  ZTensor<4> prop(2,2,2,2);
  prop = Complx(0.0, 0.0);
  double Hz = h * cos(f * t);
  double Ez = Jz / 4;
  double E  = J  / 2;
  Complx Xz(cos(Ez * dt),-sin(Ez * dt));
  Complx Ha(cos(Hz * dt),-sin(Hz * dt));
  Complx Hb(cos(Hz * dt), sin(Hz * dt));
  Complx Yz(cos(Ez * dt), sin(Ez * dt));
  Complx X (cos(E  * dt), 0.0         );
  Complx Y (0.0         ,-sin(E  * dt));
  prop(0,0,0,0) = Xz * Ha;
  prop(0,1,0,1) = Yz * X;
  prop(0,1,1,0) = Yz * Y;
  prop(1,0,0,1) = Yz * Y;
  prop(1,0,1,0) = Yz * X;
  prop(1,1,1,1) = Xz * Hb;
  return prop;
}

int real_tdmrg(ostream& fout, vector<ZTensor<3> >& mpsts, bool restart,
               int L, int M, double J, double Jz, double h, double f, int N, double dt)
{
  fout.setf(ios::fixed,ios::floatfield);
  fout.precision(4);
  fout << "====================  LATTICE INFO  ====================" << endl;
  fout << "\tL =" << setw(4) << L << ", M ="  << setw(4) << M
       << ", J =" << setw(7) << J << ", Jz =" << setw(7) << Jz << endl;

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
      if(i & 1) mpsts[i](0,1,0) = Complx( 1.0, 0.0);
      else      mpsts[i](0,0,0) = Complx( 1.0, 0.0);
    }
    mpsts[L-1].resize(M0,2,1);
    if(L & 1) mpsts[L-1](0,0,0) = Complx( 1.0, 0.0);
    else      mpsts[L-1](0,1,0) = Complx( 1.0, 0.0);
  }
  // mpo's
  vector< ZTensor<4> > mpops;
  mpops.resize(L);
  {
    // [ 0  JS+  JS-  JzSz  I ]
    ZTensor<4> mpo0(1,5,2,2);
    mpo0 = Complx( 0.0, 0.0);
    mpo0(0,1,0,1) = Complx( J /2, 0.0); // S+
    mpo0(0,2,1,0) = Complx( J /2, 0.0); // S-
    mpo0(0,3,0,0) = Complx( Jz/2, 0.0); // Sz(a)
    mpo0(0,3,1,1) = Complx(-Jz/2, 0.0); // Sz(b)
    mpo0(0,4,0,0) = Complx( 1.00, 0.0); // I(a)
    mpo0(0,4,1,1) = Complx( 1.00, 0.0); // I(b)
    Zcopy(mpo0,mpops[0]);
  }
  for(int i = 1; i < L-1; ++i) {
    // [ I   0    0    0    0 ]
    // [ S-  0    0    0    0 ]
    // [ S+  0    0    0    0 ]
    // [ Sz  0    0    0    0 ]
    // [ 0  JS+  JS-  JzSz  I ]
    ZTensor<4> mpoI(5,5,2,2);
    mpoI = Complx( 0.0, 0.0);
    mpoI(0,0,0,0) = Complx( 1.00, 0.0); // I(a)
    mpoI(0,0,1,1) = Complx( 1.00, 0.0); // I(b)
    mpoI(1,0,1,0) = Complx( 1.00, 0.0); // S-
    mpoI(2,0,0,1) = Complx( 1.00, 0.0); // S+
    mpoI(3,0,0,0) = Complx( 0.50, 0.0); // Sz(a)
    mpoI(3,0,1,1) = Complx(-0.50, 0.0); // Sz(b)
    mpoI(4,1,0,1) = Complx( J /2, 0.0); // S+
    mpoI(4,2,1,0) = Complx( J /2, 0.0); // S-
    mpoI(4,3,0,0) = Complx( Jz/2, 0.0); // Sz(a)
    mpoI(4,3,1,1) = Complx(-Jz/2, 0.0); // Sz(b)
    mpoI(4,4,0,0) = Complx( 1.00, 0.0); // I(a)
    mpoI(4,4,1,1) = Complx( 1.00, 0.0); // I(b)
    Zcopy(mpoI,mpops[i]);
  }
  {
    // [ I  ]
    // [ S- ]
    // [ S+ ]
    // [ Sz ]
    // [ 0  ]
    ZTensor<4> mpoN(5,1,2,2);
    mpoN = Complx( 0.0, 0.0);
    mpoN(0,0,0,0) = Complx( 1.00, 0.0); // I(a)
    mpoN(0,0,1,1) = Complx( 1.00, 0.0); // I(b)
    mpoN(1,0,1,0) = Complx( 1.00, 0.0); // S-
    mpoN(2,0,0,1) = Complx( 1.00, 0.0); // S+
    mpoN(3,0,0,0) = Complx( 0.50, 0.0); // Sz(a)
    mpoN(3,0,1,1) = Complx(-0.50, 0.0); // Sz(b)
    Zcopy(mpoN,mpops[L-1]);
  }

  int nodd  = L/2;
  int neven = L/2; if(L%2 == 0) --neven;
  // time-evolution
  fout << "==================== TIME-EVOLUTION ====================" << endl;
  double e   = 0.0;
  double t   = 0.0;
  int    itr = 0;
  bool   fwd = true;

//
// free-field propagator : exp(-iHt) = exp(-iEt)
//
// ZTensor<4> hodd (propagator(J, Jz, 0.0, 0.0, 0.0, dt/2));
// ZTensor<4> heven(propagator(J, Jz, 0.0, 0.0, 0.0, dt  ));
//

  while(itr < N) {
    if(itr % 100 == 0) {
      e = compute_energy(fwd,mpsts,mpops);
      fout.precision(4);
      fout << "\tt = " << setw(8) << t;
      fout.precision(16);
      fout << " : E = " << setw(20) << e << endl;
    }
    // time-evolution
    ZTensor<4> hodd (propagator(J, Jz, h, f, t, dt/2));
    ZTensor<4> heven(propagator(J, Jz, h, f, t, dt  ));
    tdmrg_sweeping(fwd,hodd,heven,mpsts,M);
    t += dt; ++itr;
  }
  e = compute_energy(fwd,mpsts,mpops);
  fout.precision(4);
  fout << "\tt = " << setw(8) << t;
  fout.precision(16);
  fout << " : E = " << setw(20) << e << endl;
  fout << "==================== FINISHED STEPS ====================" << endl;

  return 0;
}
