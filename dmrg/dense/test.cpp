#include <iostream>
#include <iomanip>

#include <vector>

#include <legacy/DENSE/TArray.h>

int main ()
{
  size_t N; std::cin >> N;

  double t = 1.0;
  double U = 1.0;
  double F = 0.5;

  std::cout.precision(8);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);

  btas::TArray<double,2> h(N,N); h.fill(0.0);

  for(size_t i = 0; i < N; ++i) {
    double v = t;
    for(size_t j = i+1; j < N; ++j) {
      h(i,j) = v;
      h(j,i) = v;
      v *= F;
    }
  }

  std::cout << "1p operator matrix :: " << h << std::endl;

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

  std::cout << "summary of SVD results :: " << q << std::endl;

  size_t K = N/2;

  std::cout << "Site[ 0] :: " << std::endl;
  std::cout << "\tui = " << std::setw(12) << q(1,0) << std::endl;

  for(size_t i = 1; i < K; ++i) {
    std::cout << "Site[" << std::setw(2) << i << "] :: " << std::endl;
    double tl = 0.0;
    double ex = 0.0;
    for(size_t j = 0; j < i; ++j) {
      tl += q(i,j)*h(i,j);
      ex += q(i,j)*q(i+1,j);
    }
    std::cout << "\ttl = " << std::setw(12) << tl << std::endl;
    std::cout << "\tex = " << std::setw(12) << ex << std::endl;
    std::cout << "\tui = " << std::setw(12) << q(i+1,i) << std::endl;
  }

  {
    size_t i = K;

    std::cout << "Site[" << std::setw(2) << i << "] :: " << std::endl;
    double tl = 0.0;
    for(size_t j = 0; j < i; ++j) tl += q(i,j)*h(i,j);
    double tr = 0.0;
    for(size_t j = i+1; j < N; ++j) tr += q(i,j)*h(i,j);

    std::cout << "\ttl = " << std::setw(12) << tl << std::endl;
    std::cout << "\tex = " << std::setw(12) << q(i,i) << std::endl;
    std::cout << "\ttr = " << std::setw(12) << tr << std::endl;
  }

  for(size_t i = K+1; i < N-1; ++i) {
    std::cout << "Site[" << std::setw(2) << i << "] :: " << std::endl;
    double ex = 0.0;
    double tr = 0.0;
    for(size_t j = i+1; j < N; ++j) {
      ex += q(i,j)*q(i-1,j);
      tr += q(i,j)*h(i,j);
    }
    std::cout << "\tvi = " << std::setw(12) << q(i-1,i) << std::endl;
    std::cout << "\tex = " << std::setw(12) << ex << std::endl;
    std::cout << "\ttr = " << std::setw(12) << tr << std::endl;
  }

  std::cout << "Site[" << std::setw(2) << N-1 << "] :: " << std::endl;
  std::cout << "\tvi = " << std::setw(12) << q(N-2,N-1) << std::endl;

  return 0;
}
