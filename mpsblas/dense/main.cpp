#include <iostream>
#include <iomanip>

#include <legacy/DENSE/TArray.h>

#include "MPX.h"
#include "hubbard.h"
#include "optimize.h"

int main ()
{
  using std::cin;
  using std::cout;
  using std::endl;
  using std::setw;

  size_t N; cin >> N;

  double t = 1.0;
  double U = 1.0;
  double F = 0.5;
  size_t M = 20;

  cout.precision(8);
  cout.setf(std::ios::fixed,std::ios::floatfield);

  btas::TArray<double,2> h(N,N); h.fill(0.0);

  for(size_t i = 0; i < N; ++i) {
    double v = t;
    for(size_t j = i+1; j < N; ++j) {
      h(i,j) = v;
      h(j,i) = v;
      v *= F;
    }
  }

  cout << "1p operator matrix :: " << h << endl;

  // Matrix q :: l = k-1, x = l-1
  // [ 0    0   0  ... ...  0  ]
  // [ u10 s11 v12 ... ... v1k ]
  // [ u20 u21 s22 v23 ... v2k ]
  // [ ... ... ... ... ... ... ]
  // [ ul0 ... ... ulx sll vlk ]
  // [  0  ... ...  0   0   0  ]

  btas::TArray<double,2> q(N,N); q.fill(0.0);

  for(size_t k = 1; k < N-1; ++k) {
    btas::TArray<double,2> W(k,N-1-k); W.fill(0.0);
    for(size_t i = 0; i < k; ++i)
      for(size_t j = k+1; j < N; ++j) W(i,j-1-k) = h(i,j);

    btas::TArray<double,1> s;
    btas::TArray<double,2> u;
    btas::TArray<double,2> v;

    btas::Gesvd('S','S',W,s,u,v);

    q(k,k) = s(0);
    for(size_t i = 0;   i < k; ++i) q(k,i) = u(i,0);
    for(size_t i = k+1; i < N; ++i) q(k,i) = v(0,i-1-k);
  }

  cout << "Summary of SVDs :: " << q << endl;

  mpsxx::MPOs<double> mpos;
  mpsxx::MPSs<double> mpss;

//mpsxx::gen_hubbard_mpos(N,mpos,t,U);

  cout << "=========== single exponential MPOs ===========" << endl;
  mpsxx::gen_hubbard_mpos(N,mpos,t,U,F);
  mpsxx::dmrg(mpos,mpss,M);

  cout << "========== LR-MPOs from Steve's idea ==========" << endl;
  mpsxx::gen_hubbard_mpos(N,mpos,h,q,U);
  mpsxx::dmrg(mpos,mpss,M);

  return 0;
}
