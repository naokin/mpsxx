#include <iostream>

#include "hubbard.h"

/// construct MPOs of Hubbard model
void mpsxx::gen_hubbard_mpos (size_t N, mpsxx::MPOs<double>& mpos, double t, double U)
{
  std::cout.precision(4);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout << "\tNearest neighbour Hubbard model :: t = " << t << ", U = " << U << std::endl;

  std::vector<size_t> D(N+1,6);
  D[0] = 1;
  D[N] = 1;

  mpos.resize(N);

  for(size_t i = 0; i < N; ++i) {
    mpos[i].resize(D[i],4,4,D[i+1]);
    mpos[i].fill(0.0);
  }

  // [ U  t*a+ -t*a-  t*b+ -t*b-  I ]

  mpos[0](0,1,1,0) =-U/2;   //  U/2 half-filling
  mpos[0](0,2,2,0) =-U/2;   //  U/2 half-filling
  mpos[0](0,3,3,0) = U;     //  U

  mpos[0](0,1,0,1) = t;     //  t a+
  mpos[0](0,3,2,1) =-t;     //  t a+ [x(-1) P(a,b)]

  mpos[0](0,0,1,2) =-t;     // -t a-
  mpos[0](0,2,3,2) = t;     // -t a- [x(-1) P(a,b)]

  mpos[0](0,2,0,3) = t;     //  t b+
  mpos[0](0,3,1,3) = t;     //  t b+

  mpos[0](0,0,2,4) =-t;     // -t b-
  mpos[0](0,1,3,4) =-t;     // -t b-

  mpos[0](0,0,0,5) = 1.0;   //  I
  mpos[0](0,1,1,5) = 1.0;   //  I
  mpos[0](0,2,2,5) = 1.0;   //  I
  mpos[0](0,3,3,5) = 1.0;   //  I

  // [ I     0      0      0      0   0 ]
  // [ a-    0      0      0      0   0 ]
  // [ a+    0      0      0      0   0 ]
  // [ b-    0      0      0      0   0 ]
  // [ b+    0      0      0      0   0 ]
  // [ U   t*a+  -t*a-   t*b+  -t*b-  I ]

  for(int i = 1; i < N-1; ++i) {

    mpos[i](0,0,0,0) = 1.0; //  I
    mpos[i](0,1,1,0) = 1.0; //  I
    mpos[i](0,2,2,0) = 1.0; //  I
    mpos[i](0,3,3,0) = 1.0; //  I

    mpos[i](1,0,1,0) = 1.0; //  a-
    mpos[i](1,2,3,0) = 1.0; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

    mpos[i](2,1,0,0) =-1.0; //  a+ [x(-1) P(l,n)]
    mpos[i](2,3,2,0) =-1.0; //  a+ [x(-1) P(a,b)]

    mpos[i](3,0,2,0) = 1.0; //  b-
    mpos[i](3,1,3,0) =-1.0; //  b- [x(-1) P(l,n)]

    mpos[i](4,2,0,0) =-1.0; //  b+ [x(-1) P(l,n)]
    mpos[i](4,3,1,0) = 1.0; //  b+

    mpos[i](5,1,1,0) =-U/2; //  U/2 half-filling
    mpos[i](5,2,2,0) =-U/2; //  U/2 half-filling
    mpos[i](5,3,3,0) = U;   //  U

    mpos[i](5,1,0,1) = t;   //  t a+
    mpos[i](5,3,2,1) =-t;   //  t a+ [x(-1) P(a,b)]

    mpos[i](5,0,1,2) =-t;   // -t a-
    mpos[i](5,2,3,2) = t;   // -t a- [x(-1) P(a,b)]

    mpos[i](5,2,0,3) = t;   //  t b+
    mpos[i](5,3,1,3) = t;   //  t b+

    mpos[i](5,0,2,4) =-t;   // -t b-
    mpos[i](5,1,3,4) =-t;   // -t b-

    mpos[i](5,0,0,5) = 1.0; //  I
    mpos[i](5,1,1,5) = 1.0; //  I
    mpos[i](5,2,2,5) = 1.0; //  I
    mpos[i](5,3,3,5) = 1.0; //  I

  }

  // [ I  ]
  // [ a- ]
  // [ a+ ]
  // [ b- ]
  // [ b+ ]
  // [ U  ]

  mpos[N-1](0,0,0,0) = 1.0; //  I
  mpos[N-1](0,1,1,0) = 1.0; //  I
  mpos[N-1](0,2,2,0) = 1.0; //  I
  mpos[N-1](0,3,3,0) = 1.0; //  I

  mpos[N-1](1,0,1,0) = 1.0; //  a-
  mpos[N-1](1,2,3,0) = 1.0; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

  mpos[N-1](2,1,0,0) =-1.0; //  a+ [x(-1) P(l,n)]
  mpos[N-1](2,3,2,0) =-1.0; //  a+ [x(-1) P(a,b)]

  mpos[N-1](3,0,2,0) = 1.0; //  b-
  mpos[N-1](3,1,3,0) =-1.0; //  b- [x(-1) P(l,n)]

  mpos[N-1](4,2,0,0) =-1.0; //  b+ [x(-1) P(l,n)]
  mpos[N-1](4,3,1,0) = 1.0; //  b+

  mpos[N-1](5,1,1,0) =-U/2; //  U/2 half-filling
  mpos[N-1](5,2,2,0) =-U/2; //  U/2 half-filling
  mpos[N-1](5,3,3,0) = U;   //  U

}

/// construct MPOs of Hubbard model having non-local single exponetial decaying terms
void mpsxx::gen_hubbard_mpos (size_t N, mpsxx::MPOs<double>& mpos, double t, double U, double F)
{
  std::vector<size_t> D(N+1,6);
  D[0] = 1;
  D[N] = 1;

  mpos.resize(N);

  for(size_t i = 0; i < N; ++i) {
    mpos[i].resize(D[i],4,4,D[i+1]);
    mpos[i].fill(0.0);
  }

  // [ U  t*a+ -t*a-  t*b+ -t*b-  I ]

  mpos[0](0,1,1,0) =-U/2;   //  U/2 half-filling
  mpos[0](0,2,2,0) =-U/2;   //  U/2 half-filling
  mpos[0](0,3,3,0) = U;     //  U

  mpos[0](0,1,0,1) = t;     //  t a+
  mpos[0](0,3,2,1) =-t;     //  t a+ [x(-1) P(a,b)]

  mpos[0](0,0,1,2) =-t;     // -t a-
  mpos[0](0,2,3,2) = t;     // -t a- [x(-1) P(a,b)]

  mpos[0](0,2,0,3) = t;     //  t b+
  mpos[0](0,3,1,3) = t;     //  t b+

  mpos[0](0,0,2,4) =-t;     // -t b-
  mpos[0](0,1,3,4) =-t;     // -t b-

  mpos[0](0,0,0,5) = 1.0;   //  I
  mpos[0](0,1,1,5) = 1.0;   //  I
  mpos[0](0,2,2,5) = 1.0;   //  I
  mpos[0](0,3,3,5) = 1.0;   //  I

  // [ I     0      0      0      0   0 ]
  // [ a-  F*I      0      0      0   0 ]
  // [ a+    0    F*I      0      0   0 ]
  // [ b-    0      0    F*I      0   0 ]
  // [ b+    0      0      0    F*I   0 ]
  // [ U   t*a+  -t*a-   t*b+  -t*b-  I ]

  for(int i = 1; i < N-1; ++i) {

    mpos[i](0,0,0,0) = 1.0; //  I
    mpos[i](0,1,1,0) = 1.0; //  I
    mpos[i](0,2,2,0) = 1.0; //  I
    mpos[i](0,3,3,0) = 1.0; //  I

    mpos[i](1,0,1,0) = 1.0; //  a-
    mpos[i](1,2,3,0) = 1.0; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

    mpos[i](1,0,0,1) = F;   //  F I
    mpos[i](1,1,1,1) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](1,2,2,1) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](1,3,3,1) = F;   //  F I

    mpos[i](2,1,0,0) =-1.0; //  a+ [x(-1) P(l,n)]
    mpos[i](2,3,2,0) =-1.0; //  a+ [x(-1) P(a,b)]

    mpos[i](2,0,0,2) = F;   //  F I
    mpos[i](2,1,1,2) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](2,2,2,2) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](2,3,3,2) = F;   //  F I

    mpos[i](3,0,2,0) = 1.0; //  b-
    mpos[i](3,1,3,0) =-1.0; //  b- [x(-1) P(l,n)]

    mpos[i](3,0,0,3) = F;   //  F I
    mpos[i](3,1,1,3) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](3,2,2,3) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](3,3,3,3) = F;   //  F I

    mpos[i](4,2,0,0) =-1.0; //  b+ [x(-1) P(l,n)]
    mpos[i](4,3,1,0) = 1.0; //  b+

    mpos[i](4,0,0,4) = F;   //  F I
    mpos[i](4,1,1,4) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](4,2,2,4) =-F;   //  F I [x(-1) P(l,n)]
    mpos[i](4,3,3,4) = F;   //  F I

    mpos[i](5,1,1,0) =-U/2; //  U/2 half-filling
    mpos[i](5,2,2,0) =-U/2; //  U/2 half-filling
    mpos[i](5,3,3,0) = U;   //  U

    mpos[i](5,1,0,1) = t;   //  t a+
    mpos[i](5,3,2,1) =-t;   //  t a+ [x(-1) P(a,b)]

    mpos[i](5,0,1,2) =-t;   // -t a-
    mpos[i](5,2,3,2) = t;   // -t a- [x(-1) P(a,b)]

    mpos[i](5,2,0,3) = t;   //  t b+
    mpos[i](5,3,1,3) = t;   //  t b+

    mpos[i](5,0,2,4) =-t;   // -t b-
    mpos[i](5,1,3,4) =-t;   // -t b-

    mpos[i](5,0,0,5) = 1.0; //  I
    mpos[i](5,1,1,5) = 1.0; //  I
    mpos[i](5,2,2,5) = 1.0; //  I
    mpos[i](5,3,3,5) = 1.0; //  I

  }

  // [ I  ]
  // [ a- ]
  // [ a+ ]
  // [ b- ]
  // [ b+ ]
  // [ U  ]

  mpos[N-1](0,0,0,0) = 1.0; //  I
  mpos[N-1](0,1,1,0) = 1.0; //  I
  mpos[N-1](0,2,2,0) = 1.0; //  I
  mpos[N-1](0,3,3,0) = 1.0; //  I

  mpos[N-1](1,0,1,0) = 1.0; //  a-
  mpos[N-1](1,2,3,0) = 1.0; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

  mpos[N-1](2,1,0,0) =-1.0; //  a+ [x(-1) P(l,n)]
  mpos[N-1](2,3,2,0) =-1.0; //  a+ [x(-1) P(a,b)]

  mpos[N-1](3,0,2,0) = 1.0; //  b-
  mpos[N-1](3,1,3,0) =-1.0; //  b- [x(-1) P(l,n)]

  mpos[N-1](4,2,0,0) =-1.0; //  b+ [x(-1) P(l,n)]
  mpos[N-1](4,3,1,0) = 1.0; //  b+

  mpos[N-1](5,1,1,0) =-U/2; //  U/2 half-filling
  mpos[N-1](5,2,2,0) =-U/2; //  U/2 half-filling
  mpos[N-1](5,3,3,0) = U;   //  U

}

/// construct MPOs of Hubbard model having non-local single exponetial decaying terms (Steve's idea)
void mpsxx::gen_hubbard_mpos (size_t N, mpsxx::MPOs<double>& mpos, const btas::TArray<double,2>& t, const btas::TArray<double,2>& q, double U)
{
  std::vector<size_t> D(N+1,6);
  D[0] = 1;
  D[N] = 1;

  mpos.resize(N);

  for(size_t i = 0; i < N; ++i) {
    mpos[i].resize(D[i],4,4,D[i+1]);
    mpos[i].fill(0.0);
  }

  // the first site

  // [ U  t*a+ -t*a-  t*b+ -t*b-  I ]

  {
    double ui = q(1,0);

    mpos[0](0,1,1,0) =-U/2;   //  U/2 half-filling
    mpos[0](0,2,2,0) =-U/2;   //  U/2 half-filling
    mpos[0](0,3,3,0) = U;     //  U

    mpos[0](0,1,0,1) = ui;     //  t a+
    mpos[0](0,3,2,1) =-ui;     //  t a+ [x(-1) P(a,b)]

    mpos[0](0,0,1,2) =-ui;     // -t a-
    mpos[0](0,2,3,2) = ui;     // -t a- [x(-1) P(a,b)]

    mpos[0](0,2,0,3) = ui;     //  t b+
    mpos[0](0,3,1,3) = ui;     //  t b+

    mpos[0](0,0,2,4) =-ui;     // -t b-
    mpos[0](0,1,3,4) =-ui;     // -t b-

    mpos[0](0,0,0,5) = 1.0;   //  I
    mpos[0](0,1,1,5) = 1.0;   //  I
    mpos[0](0,2,2,5) = 1.0;   //  I
    mpos[0](0,3,3,5) = 1.0;   //  I
  }

  // [ I     0      0      0      0   0 ]
  // [ a-  F*I      0      0      0   0 ]
  // [ a+    0    F*I      0      0   0 ]
  // [ b-    0      0    F*I      0   0 ]
  // [ b+    0      0      0    F*I   0 ]
  // [ U   t*a+  -t*a-   t*b+  -t*b-  I ]

  size_t K = N/2;

  // Left sites

  for(int i = 1; i < K; ++i) {

    double tl = 0.0;
    double ex = 0.0;
    for(int j = 0; j < i; ++j) {
      tl += q(i,j)*t(i,j);
      ex += q(i,j)*q(i+1,j);
    }
    double ui = q(i+1,i);

    mpos[i](0,0,0,0) = 1.0; //  I
    mpos[i](0,1,1,0) = 1.0; //  I
    mpos[i](0,2,2,0) = 1.0; //  I
    mpos[i](0,3,3,0) = 1.0; //  I

    mpos[i](1,0,1,0) = tl; //  a-
    mpos[i](1,2,3,0) = tl; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

    mpos[i](1,0,0,1) = ex;   //  F I
    mpos[i](1,1,1,1) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](1,2,2,1) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](1,3,3,1) = ex;   //  F I

    mpos[i](2,1,0,0) =-tl; //  a+ [x(-1) P(l,n)]
    mpos[i](2,3,2,0) =-tl; //  a+ [x(-1) P(a,b)]

    mpos[i](2,0,0,2) = ex;   //  F I
    mpos[i](2,1,1,2) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](2,2,2,2) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](2,3,3,2) = ex;   //  F I

    mpos[i](3,0,2,0) = tl; //  b-
    mpos[i](3,1,3,0) =-tl; //  b- [x(-1) P(l,n)]

    mpos[i](3,0,0,3) = ex;   //  F I
    mpos[i](3,1,1,3) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](3,2,2,3) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](3,3,3,3) = ex;   //  F I

    mpos[i](4,2,0,0) =-tl; //  b+ [x(-1) P(l,n)]
    mpos[i](4,3,1,0) = tl; //  b+

    mpos[i](4,0,0,4) = ex;   //  F I
    mpos[i](4,1,1,4) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](4,2,2,4) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](4,3,3,4) = ex;   //  F I

    mpos[i](5,1,1,0) =-U/2; //  U/2 half-filling
    mpos[i](5,2,2,0) =-U/2; //  U/2 half-filling
    mpos[i](5,3,3,0) = U;   //  U

    mpos[i](5,1,0,1) = ui;   //  t a+
    mpos[i](5,3,2,1) =-ui;   //  t a+ [x(-1) P(a,b)]

    mpos[i](5,0,1,2) =-ui;   // -t a-
    mpos[i](5,2,3,2) = ui;   // -t a- [x(-1) P(a,b)]

    mpos[i](5,2,0,3) = ui;   //  t b+
    mpos[i](5,3,1,3) = ui;   //  t b+

    mpos[i](5,0,2,4) =-ui;   // -t b-
    mpos[i](5,1,3,4) =-ui;   // -t b-

    mpos[i](5,0,0,5) = 1.0; //  I
    mpos[i](5,1,1,5) = 1.0; //  I
    mpos[i](5,2,2,5) = 1.0; //  I
    mpos[i](5,3,3,5) = 1.0; //  I

  }

  // Central sites

  {

    int i = K;

    double tl = 0.0;
    for(int j = 0; j < i; ++j) {
      tl += q(i,j)*t(i,j);
    }
    double tr = 0.0;
    for(int j = i+1; j < N; ++j) {
      tr += q(i,j)*t(i,j);
    }
    double ex = q(i,i);

    mpos[i](0,0,0,0) = 1.0; //  I
    mpos[i](0,1,1,0) = 1.0; //  I
    mpos[i](0,2,2,0) = 1.0; //  I
    mpos[i](0,3,3,0) = 1.0; //  I

    mpos[i](1,0,1,0) = tl; //  a-
    mpos[i](1,2,3,0) = tl; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

    mpos[i](1,0,0,1) = ex;   //  F I
    mpos[i](1,1,1,1) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](1,2,2,1) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](1,3,3,1) = ex;   //  F I

    mpos[i](2,1,0,0) =-tl; //  a+ [x(-1) P(l,n)]
    mpos[i](2,3,2,0) =-tl; //  a+ [x(-1) P(a,b)]

    mpos[i](2,0,0,2) = ex;   //  F I
    mpos[i](2,1,1,2) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](2,2,2,2) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](2,3,3,2) = ex;   //  F I

    mpos[i](3,0,2,0) = tl; //  b-
    mpos[i](3,1,3,0) =-tl; //  b- [x(-1) P(l,n)]

    mpos[i](3,0,0,3) = ex;   //  F I
    mpos[i](3,1,1,3) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](3,2,2,3) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](3,3,3,3) = ex;   //  F I

    mpos[i](4,2,0,0) =-tl; //  b+ [x(-1) P(l,n)]
    mpos[i](4,3,1,0) = tl; //  b+

    mpos[i](4,0,0,4) = ex;   //  F I
    mpos[i](4,1,1,4) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](4,2,2,4) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](4,3,3,4) = ex;   //  F I

    mpos[i](5,1,1,0) =-U/2; //  U/2 half-filling
    mpos[i](5,2,2,0) =-U/2; //  U/2 half-filling
    mpos[i](5,3,3,0) = U;   //  U

    mpos[i](5,1,0,1) = tr;   //  t a+
    mpos[i](5,3,2,1) =-tr;   //  t a+ [x(-1) P(a,b)]

    mpos[i](5,0,1,2) =-tr;   // -t a-
    mpos[i](5,2,3,2) = tr;   // -t a- [x(-1) P(a,b)]

    mpos[i](5,2,0,3) = tr;   //  t b+
    mpos[i](5,3,1,3) = tr;   //  t b+

    mpos[i](5,0,2,4) =-tr;   // -t b-
    mpos[i](5,1,3,4) =-tr;   // -t b-

    mpos[i](5,0,0,5) = 1.0; //  I
    mpos[i](5,1,1,5) = 1.0; //  I
    mpos[i](5,2,2,5) = 1.0; //  I
    mpos[i](5,3,3,5) = 1.0; //  I

  }

  // Right sites

  for(int i = K+1; i < N-1; ++i) {

    double tr = 0.0;
    double ex = 0.0;
    for(int j = i+1; j < N; ++j) {
      tr += q(i,j)*t(i,j);
      ex += q(i,j)*q(i-1,j);
    }
    double vi = q(i-1,i);

    mpos[i](0,0,0,0) = 1.0; //  I
    mpos[i](0,1,1,0) = 1.0; //  I
    mpos[i](0,2,2,0) = 1.0; //  I
    mpos[i](0,3,3,0) = 1.0; //  I

    mpos[i](1,0,1,0) = vi; //  a-
    mpos[i](1,2,3,0) = vi; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

    mpos[i](1,0,0,1) = ex;   //  F I
    mpos[i](1,1,1,1) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](1,2,2,1) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](1,3,3,1) = ex;   //  F I

    mpos[i](2,1,0,0) =-vi; //  a+ [x(-1) P(l,n)]
    mpos[i](2,3,2,0) =-vi; //  a+ [x(-1) P(a,b)]

    mpos[i](2,0,0,2) = ex;   //  F I
    mpos[i](2,1,1,2) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](2,2,2,2) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](2,3,3,2) = ex;   //  F I

    mpos[i](3,0,2,0) = vi; //  b-
    mpos[i](3,1,3,0) =-vi; //  b- [x(-1) P(l,n)]

    mpos[i](3,0,0,3) = ex;   //  F I
    mpos[i](3,1,1,3) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](3,2,2,3) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](3,3,3,3) = ex;   //  F I

    mpos[i](4,2,0,0) =-vi; //  b+ [x(-1) P(l,n)]
    mpos[i](4,3,1,0) = vi; //  b+

    mpos[i](4,0,0,4) = ex;   //  F I
    mpos[i](4,1,1,4) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](4,2,2,4) =-ex;   //  F I [x(-1) P(l,n)]
    mpos[i](4,3,3,4) = ex;   //  F I

    mpos[i](5,1,1,0) =-U/2; //  U/2 half-filling
    mpos[i](5,2,2,0) =-U/2; //  U/2 half-filling
    mpos[i](5,3,3,0) = U;   //  U

    mpos[i](5,1,0,1) = tr;   //  t a+
    mpos[i](5,3,2,1) =-tr;   //  t a+ [x(-1) P(a,b)]

    mpos[i](5,0,1,2) =-tr;   // -t a-
    mpos[i](5,2,3,2) = tr;   // -t a- [x(-1) P(a,b)]

    mpos[i](5,2,0,3) = tr;   //  t b+
    mpos[i](5,3,1,3) = tr;   //  t b+

    mpos[i](5,0,2,4) =-tr;   // -t b-
    mpos[i](5,1,3,4) =-tr;   // -t b-

    mpos[i](5,0,0,5) = 1.0; //  I
    mpos[i](5,1,1,5) = 1.0; //  I
    mpos[i](5,2,2,5) = 1.0; //  I
    mpos[i](5,3,3,5) = 1.0; //  I

  }

  // [ I  ]
  // [ a- ]
  // [ a+ ]
  // [ b- ]
  // [ b+ ]
  // [ U  ]

  {
    double vi = q(N-2,N-1);

    mpos[N-1](0,0,0,0) = 1.0; //  I
    mpos[N-1](0,1,1,0) = 1.0; //  I
    mpos[N-1](0,2,2,0) = 1.0; //  I
    mpos[N-1](0,3,3,0) = 1.0; //  I

    mpos[N-1](1,0,1,0) = vi; //  a-
    mpos[N-1](1,2,3,0) = vi; //  a- [x(-1) P(a,b)] [x(-1) P(l,n)]

    mpos[N-1](2,1,0,0) =-vi; //  a+ [x(-1) P(l,n)]
    mpos[N-1](2,3,2,0) =-vi; //  a+ [x(-1) P(a,b)]

    mpos[N-1](3,0,2,0) = vi; //  b-
    mpos[N-1](3,1,3,0) =-vi; //  b- [x(-1) P(l,n)]

    mpos[N-1](4,2,0,0) =-vi; //  b+ [x(-1) P(l,n)]
    mpos[N-1](4,3,1,0) = vi; //  b+

    mpos[N-1](5,1,1,0) =-U/2; //  U/2 half-filling
    mpos[N-1](5,2,2,0) =-U/2; //  U/2 half-filling
    mpos[N-1](5,3,3,0) = U;   //  U
  }
}
