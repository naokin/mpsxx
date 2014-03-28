#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <dense.h>
using namespace std;
using namespace btas;

#define DEBUG(msg) { cout << "DEBUG:" << msg << endl; }

void canonicalize(const bool& forward, DTensor<3>& mps, DTensor<2>& res)
{
  if(forward) {
    DTensor<1> s;
    DTensor<3> u;
    Dgesvd(mps,shape(2),s,u,res);
    Dcopy(u,mps);
    Dright_update(s,res);
  }
  else {
    DTensor<1> s;
    DTensor<3> v;
    Dgesvd(mps,shape(1,2),s,res,v);
    Dcopy(v, mps);
    Dleft_update(res,s);
  }
}

void renormalize(const bool& forward, const DTensor<3>& mps, const DTensor<4>& mpo, const DTensor<3>& stin, DTensor<3>& stout)
{
  stout.free();
  if(forward) {
    DTensor<4> scr1;
    Dcontract(1.0,stin,shape(0),mps,shape(0),1.0,scr1);
    DTensor<4> scr2;
    Dcontract(1.0,scr1,shape(0,2),mpo,shape(0,2),1.0,scr2);
    Dcontract(1.0,scr2,shape(0,3),mps,shape(0,1),1.0,stout);
  }
  else {
    DTensor<4> scr1;
    Dcontract(1.0,mps,shape(2),stin,shape(0),1.0,scr1);
    DTensor<4> scr2;
    Dcontract(1.0,scr1,shape(1,2),mpo,shape(2,1),1.0,scr2);
    Dcontract(1.0,scr2,shape(3,1),mps,shape(1,2),1.0,stout);
  }
}

void h_diagonal(const DTensor<4>& mpo, const DTensor<3>& left, const DTensor<3>& right, DTensor<3>& hdiag)
{
  int ds = mpo.extent(2);
  int Dl = mpo.extent(0);
  int Dr = mpo.extent(1);
  int Ml = left.extent(0);
  int Mr = right.extent(0);
  DTensor<3> dmpo(Dl,Dr,ds);
  for(int i = 0; i < Dl; ++i) {
    for(int j = 0; j < Dr; ++j) {
      for(int k = 0; k < ds; ++k) {
        dmpo(i,j,k) = mpo(i,j,k,k);
      }
    }
  }
  DTensor<2> dleft(Ml,Dl);
  for(int i = 0; i < Ml; ++i) {
    for(int j = 0; j < Dl; ++j) {
      dleft(i,j) = left(i,j,i);
    }
  }
  DTensor<2> dright(Mr,Dr);
  for(int i = 0; i < Mr; ++i) {
    for(int j = 0; j < Dr; ++j) {
      dright(i,j) = right(i,j,i);
    }
  }
  DTensor<3> scr;
  Dcontract(1.0,dleft,shape(1),dmpo,shape(0),1.0,scr);
  hdiag.free();
  Dcontract(1.0,scr,shape(1),dright,shape(1),1.0,hdiag);
}

void sigmavector(const DTensor<3>& mps, const DTensor<4>& mpo, const DTensor<3>& left, const DTensor<3>& right, DTensor<3>& sgv)
{
  DTensor<4> scr1;
  Dcontract(1.0,left,shape(2),mps,shape(0),1.0,scr1);
  DTensor<4> scr2;
  Dcontract(1.0,scr1,shape(1,2),mpo,shape(0,3),1.0,scr2);
  sgv.free();
  Dcontract(1.0,scr2,shape(2,1),right,shape(1,2),1.0,sgv);
}

template<int N>
void normalize(DTensor<N>& t)
{
  double norm = Ddot(t, t);
  Dscal(1.0/sqrt(norm),t);
}

template<int N>
void orthogonalize(const DTensor<N>& base, DTensor<N>& sub)
{
  double ovlp = Ddot(base, sub);
  Daxpy(-ovlp, base, sub);
}

void precondition(const double& eig, const DTensor<3>& hdiag, DTensor<3>& mps)
{
  DTensor<3>::const_iterator it_diag = hdiag.begin();
  DTensor<3>::iterator it_mps = mps.begin();
  for(; it_mps != mps.end(); ++it_diag, ++it_mps) {
    double denom = eig - *it_diag;
    *it_mps /= denom;
  }
}

double davidson(DTensor<3>& mps, const DTensor<4>& mpo, const DTensor<3>& left, const DTensor<3>& right, const DTensor<3>& hdiag)
{
  const int MXITER = 100;
  const int MXRVEC = 20;
  vector< DTensor<3> > trial_vec(MXRVEC,DTensor<3>());
  vector< DTensor<3> > sigma_vec(MXRVEC,DTensor<3>());

  normalize(mps);
  Dcopy(mps,trial_vec[0]);

  int iter = 0;
  bool conv = false;
  double energy = 0.0;
  while(!conv && iter < MXITER) {
    DTensor<3> eigen_vec;
    DTensor<3> error_vec;
    for(int m = 1; !conv && m < MXRVEC; ++m) {
      sigmavector(trial_vec[m-1],mpo,left,right,sigma_vec[m-1]);

      DTensor<2> dav_small(m, m);
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          dav_small(i, j) = Ddot(sigma_vec[i],trial_vec[j]);
        }
      }
      DTensor<2> dav_egvec(m, m);
      DTensor<1> dav_egval(m);
      Dsyev(dav_small,dav_egval,dav_egvec);
      energy = dav_egval(0);

      vector< DTensor<3> > trial_vec_save(m,DTensor<3>());
      vector< DTensor<3> > sigma_vec_save(m,DTensor<3>());

      for(int i = 0; i < m; ++i) {
        Dcopy(trial_vec[i],trial_vec_save[i]);
        Dcopy(sigma_vec[i],sigma_vec_save[i]);
        Dscal(dav_egvec(i,i),trial_vec[i]);
        Dscal(dav_egvec(i,i),sigma_vec[i]);
      }
      for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
          if(i != j) {
            Daxpy(dav_egvec(i,j),trial_vec_save[j],trial_vec[i]);
            Daxpy(dav_egvec(i,j),sigma_vec_save[j],sigma_vec[i]);
          }
        }
      }
      Dcopy(trial_vec[0],eigen_vec);
      Dcopy(sigma_vec[0],error_vec);
      Daxpy(-energy,eigen_vec,error_vec);
      double radius = Ddot(error_vec,error_vec);
      if(radius < 1.0e-10) {
        conv = true;
        Dcopy(eigen_vec,mps);
      }
      else {
        precondition(energy,hdiag,error_vec);
        for(int i = 0; i < m; ++i) {
          normalize(error_vec);
          orthogonalize(trial_vec[i],error_vec);
        }
        normalize(error_vec);
        Dcopy(error_vec,trial_vec[m]);
      }
    }
    if(!conv) {
      for(int i = 0; i < MXRVEC; ++i) {
        trial_vec[i].free();
        sigma_vec[i].free();
      }
      Dcopy(eigen_vec,trial_vec[0]);
      ++iter;
    }
  }
  return energy;
}

double optimize(DTensor<3>& mps, const DTensor<4>& mpo, const DTensor<3>& left, const DTensor<3>& right)
{
  DTensor<3> hdiag;
  h_diagonal(mpo,left,right,hdiag);
  double energy = davidson(mps,mpo,left,right,hdiag);
  return energy;
}

int dmrg(ostream& fout, vector<DTensor<3> >& mpsts, bool restart, int L, int M, double J, double Jz)
{
  fout.setf(ios::fixed,ios::floatfield);
  fout.precision(4);
  fout << "====================  LATTICE INFO  ====================" << endl;
  fout << "\tL =" << setw(4) << L << ", M ="  << setw(4) << M
       << ", J =" << setw(7) << J << ", Jz =" << setw(7) << Jz << endl;

  DTensor<4> mpo0(1,5,2,2);
  DTensor<4> mpoI(5,5,2,2);
  DTensor<4> mpoN(5,1,2,2);

  // [ 0  JS+  JS-  JzSz  I ]
  mpo0 = 0.0;
  mpo0(0,1,0,1) = J /2; // S+
  mpo0(0,2,1,0) = J /2; // S-
  mpo0(0,3,0,0) = Jz/2; // Sz(a)
  mpo0(0,3,1,1) =-Jz/2; // Sz(b)
  mpo0(0,4,0,0) = 1.00; // I(a)
  mpo0(0,4,1,1) = 1.00; // I(b)

  // [ I   0    0    0    0 ]
  // [ S-  0    0    0    0 ]
  // [ S+  0    0    0    0 ]
  // [ Sz  0    0    0    0 ]
  // [ 0  JS+  JS-  JzSz  I ]
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

  // [ I  ]
  // [ S- ]
  // [ S+ ]
  // [ Sz ]
  // [ 0  ]
  mpoN = 0.0;
  mpoN(0,0,0,0) = 1.00; // I(a)
  mpoN(0,0,1,1) = 1.00; // I(b)
  mpoN(1,0,1,0) = 1.00; // S-
  mpoN(2,0,0,1) = 1.00; // S+
  mpoN(3,0,0,0) = 0.50; // Sz(a)
  mpoN(3,0,1,1) =-0.50; // Sz(b)

  // dummy storage
  DTensor<3> stdum(1,1,1);
  stdum = 1.0;

  // initialize mpsts
  if(!restart) {
    mpsts.resize(L);
    mpsts[0].resize(1,2,M);
    for(int i = 1; i < L-1; ++i) {
      mpsts[i].resize(M,2,M);
    }
    mpsts[L-1].resize(M,2,1);

    srand(time(NULL));
    for(int i = 0; i < L; ++i) {
      for(DTensor<3>::iterator it = mpsts[i].begin(); it != mpsts[i].end(); ++it) {
        *it = 2.0 * ((double) rand() / RAND_MAX - 0.5);
      }
    }
  }
  // initialize storages
  vector< DTensor<3> > storages;
  storages.resize(L);
  // set thermodynamic limit mps and backward canonicalization
  for(int i = L-1; i > 0; --i) {
    DTensor<2> res;
    canonicalize(0,mpsts[i],res);
    DTensor<3> guess;
    Dgemm(NoTrans,NoTrans,1.0,mpsts[i-1],res,1.0,guess);
    Dcopy(guess,mpsts[i-1]);
    if(i == L-1) {
      renormalize(0,mpsts[i],mpoN,stdum,storages[i]);
    }
    else {
      renormalize(0,mpsts[i],mpoI,storages[i+1],storages[i]);
    }
  }
  normalize(mpsts[0]);

  // starting sweep optimization
  fout.precision(16);
  double energy = 0.0;
  double radius = 1.0;
  int iter = 0;
  while(iter < 100 && radius > 1.0e-8) {
    fout << "==================== SWEEP ITER[" << setw(2) << iter << "] ====================" << endl;
    double energy_save = energy;
    vector<double> sweep_energy;
    sweep_energy.reserve(2*L);
    // forward sweep
    fout << "++++++++++++++++++++ FORWARD  SWEEP ++++++++++++++++++++" << endl;
    for(int i = 0; i < L-1; ++i) {
      double esite;
      if(i == 0) {
        esite = optimize(mpsts[i],mpo0,stdum,storages[i+1]);
        DTensor<2> res;
        canonicalize(1,mpsts[i],res);
        DTensor<3> guess;
        Dgemm(NoTrans,NoTrans,1.0,res,mpsts[i+1],1.0,guess);
        Dcopy(guess,mpsts[i+1]);
        renormalize(1,mpsts[i],mpo0,stdum,storages[i]);
      }
      else {
        esite = optimize(mpsts[i],mpoI,storages[i-1],storages[i+1]);
        DTensor<2> res;
        canonicalize(1,mpsts[i],res);
        DTensor<3> guess;
        Dgemm(NoTrans,NoTrans,1.0,res,mpsts[i+1],1.0,guess);
        Dcopy(guess,mpsts[i+1]);
        renormalize(1,mpsts[i],mpoI,storages[i-1],storages[i]);
      }
      fout << "\tENERGY @ SITE[" << setw(3) << i << "] = " << setw(20) << esite << endl;
      sweep_energy.push_back(esite);
    }
    // backward sweep
    fout << "++++++++++++++++++++ BACKWARD SWEEP ++++++++++++++++++++" << endl;
    for(int i = L-1; i > 0; --i) {
      double esite;
      if(i == L-1) {
        esite = optimize(mpsts[i],mpoN,storages[i-1],stdum);
        DTensor<2> res;
        canonicalize(0,mpsts[i],res);
        DTensor<3> guess;
        Dgemm(NoTrans,NoTrans,1.0,mpsts[i-1],res,1.0,guess);
        Dcopy(guess,mpsts[i-1]);
        renormalize(0,mpsts[i],mpoN,stdum,storages[i]);
      }
      else {
        esite = optimize(mpsts[i],mpoI,storages[i-1],storages[i+1]);
        DTensor<2> res;
        canonicalize(0,mpsts[i],res);
        DTensor<3> guess;
        Dgemm(NoTrans,NoTrans,1.0,mpsts[i-1],res,1.0,guess);
        Dcopy(guess,mpsts[i-1]);
        renormalize(0,mpsts[i],mpoI,storages[i+1],storages[i]);
      }
      fout << "\tENERGY @ SITE[" << setw(3) << i << "] = " << setw(20) << esite << endl;
      sweep_energy.push_back(esite);
    }
    energy = *min_element(sweep_energy.begin(),sweep_energy.end());
    radius = fabs(energy-energy_save);
    sweep_energy.clear();
    ++iter;
  }
  fout << "==================== FINISHED SWEEP ====================" << endl;

  return 0;
}
