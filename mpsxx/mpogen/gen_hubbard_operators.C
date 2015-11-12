#include <iostream>
#include <iomanip>

#include <vector>

#include <legacy/QSPARSE/QSTArray.h>

#include <MpSite.h>
#include <driver/fileio.h>

#include "gen_hubbard_operators.h"
#include "gen_site_operator.h"
#include "prime_operators.h"

void mpsxx::gen_hubbard_operators
(mpsxx::MPOs<double, fermion>& mpos, const double& t, const double& u, const std::string& prefix)
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;

  size_t N = mpos.size();

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tCONSTRUCT FERMIONIC-HUBBARD MATRIX PRODUCT OPERATORS (MPOs) "                                       << endl;
  cout.precision(4);
  cout << "\t\t\t+ coupling coefficient t  : " << setw(8) << fixed << t  << endl; 
  cout << "\t\t\t+ coupling coefficient U  : " << setw(8) << fixed << u  << endl;
  cout << "\t====================================================================================================" << endl;

  btas::Qshapes<fermion> qz; // 0 quantum number
  qz.push_back(fermion( 0,  0));

  btas::Qshapes<fermion> qi; // quantum index comes in
  qi.push_back(fermion( 0,  0)); // 0
  qi.push_back(fermion( 1,  1)); // a-
  qi.push_back(fermion(-1, -1)); // a+
  qi.push_back(fermion( 1, -1)); // b-
  qi.push_back(fermion(-1,  1)); // b+
  qi.push_back(fermion( 0,  0)); // I

  btas::Qshapes<fermion> qo; // quantum index comes out
  qo.push_back(fermion( 0,  0)); // I
  qo.push_back(fermion(-1, -1)); // a+
  qo.push_back(fermion( 1,  1)); // a-
  qo.push_back(fermion(-1,  1)); // b+
  qo.push_back(fermion( 1, -1)); // b-
  qo.push_back(fermion( 0,  0)); // 0

  // resize & set to 0
  mpos[ 0 ].resize(fermion::zero(), make_array(qz, MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(), qo));
  for(int i = 1; i < N-1; ++i)
    mpos[i].resize(fermion::zero(), make_array(qi, MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(), qo));
  mpos[N-1].resize(fermion::zero(), make_array(qi, MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(), qz));

  // set block elements
  btas::TArray<double,4> data_Ip(1, 1, 1, 1); data_Ip = 1.0;
  btas::TArray<double,4> data_Im(1, 1, 1, 1); data_Im =-1.0;
  btas::TArray<double,4> data_tp(1, 1, 1, 1); data_tp = t;
  btas::TArray<double,4> data_tm(1, 1, 1, 1); data_tm =-t;
  btas::TArray<double,4> data_Un(1, 1, 1, 1); data_Un = u;
  // insert blocks
  mpos[ 0 ].insert(btas::shape(0, 3, 3, 0), data_Un); //  U
  mpos[ 0 ].insert(btas::shape(0, 1, 0, 1), data_tp); //  t a+
  mpos[ 0 ].insert(btas::shape(0, 3, 2, 1), data_tm); //  t a+ [x(-1) due to P(a,b)]
  mpos[ 0 ].insert(btas::shape(0, 0, 1, 2), data_tm); // -t a-
  mpos[ 0 ].insert(btas::shape(0, 2, 3, 2), data_tp); // -t a- [x(-1) due to P(a,b)]
  mpos[ 0 ].insert(btas::shape(0, 2, 0, 3), data_tp); //  t b+
  mpos[ 0 ].insert(btas::shape(0, 3, 1, 3), data_tp); //  t b+
  mpos[ 0 ].insert(btas::shape(0, 0, 2, 4), data_tm); // -t b-
  mpos[ 0 ].insert(btas::shape(0, 1, 3, 4), data_tm); // -t b-
  mpos[ 0 ].insert(btas::shape(0, 0, 0, 5), data_Ip); //  I
  mpos[ 0 ].insert(btas::shape(0, 1, 1, 5), data_Ip); //  I
  mpos[ 0 ].insert(btas::shape(0, 2, 2, 5), data_Ip); //  I
  mpos[ 0 ].insert(btas::shape(0, 3, 3, 5), data_Ip); //  I
  for(int i = 1; i < N-1; ++i) {
    mpos[i].insert(btas::shape(0, 0, 0, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(0, 1, 1, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(0, 2, 2, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(0, 3, 3, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(1, 0, 1, 0), data_Ip); //  a-
    mpos[i].insert(btas::shape(1, 2, 3, 0), data_Im); //  a- [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(2, 1, 0, 0), data_Ip); //  a+
    mpos[i].insert(btas::shape(2, 3, 2, 0), data_Im); //  a+ [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(3, 0, 2, 0), data_Ip); //  b-
    mpos[i].insert(btas::shape(3, 1, 3, 0), data_Ip); //  b-
    mpos[i].insert(btas::shape(4, 2, 0, 0), data_Ip); //  b+
    mpos[i].insert(btas::shape(4, 3, 1, 0), data_Ip); //  b+
    mpos[i].insert(btas::shape(5, 3, 3, 0), data_Un); //  U
    mpos[i].insert(btas::shape(5, 1, 0, 1), data_tp); //  t a+
    mpos[i].insert(btas::shape(5, 3, 2, 1), data_tm); //  t a+ [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(5, 0, 1, 2), data_tm); // -t a-
    mpos[i].insert(btas::shape(5, 2, 3, 2), data_tp); // -t a- [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(5, 2, 0, 3), data_tp); //  t b+
    mpos[i].insert(btas::shape(5, 3, 1, 3), data_tp); //  t b+
    mpos[i].insert(btas::shape(5, 0, 2, 4), data_tm); // -t b-
    mpos[i].insert(btas::shape(5, 1, 3, 4), data_tm); // -t b-
    mpos[i].insert(btas::shape(5, 0, 0, 5), data_Ip); //  I
    mpos[i].insert(btas::shape(5, 1, 1, 5), data_Ip); //  I
    mpos[i].insert(btas::shape(5, 2, 2, 5), data_Ip); //  I
    mpos[i].insert(btas::shape(5, 3, 3, 5), data_Ip); //  I
  }
  mpos[N-1].insert(btas::shape(0, 0, 0, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(0, 1, 1, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(0, 2, 2, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(0, 3, 3, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(1, 0, 1, 0), data_Ip); //  a-
  mpos[N-1].insert(btas::shape(1, 2, 3, 0), data_Im); //  a- [x(-1) due to P(a,b)]
  mpos[N-1].insert(btas::shape(2, 1, 0, 0), data_Ip); //  a+
  mpos[N-1].insert(btas::shape(2, 3, 2, 0), data_Im); //  a+ [x(-1) due to P(a,b)]
  mpos[N-1].insert(btas::shape(3, 0, 2, 0), data_Ip); //  b-
  mpos[N-1].insert(btas::shape(3, 1, 3, 0), data_Ip); //  b-
  mpos[N-1].insert(btas::shape(4, 2, 0, 0), data_Ip); //  b+
  mpos[N-1].insert(btas::shape(4, 3, 1, 0), data_Ip); //  b+
  mpos[N-1].insert(btas::shape(5, 3, 3, 0), data_Un); //  U

  // taking operator parity P(O(l),<n|)
  std::vector<int> indx1(1, 0); // left mpo index [0]
  std::vector<int> indx2(1, 1); // phys bra index [1]
  for(int i = 0; i < N; ++i) {
    mpos[i].check_dshape();
    mpos[i].parity(indx1, indx2);
    save(mpos[i], get_mpofile(prefix, i));
    mpos[i].clear();
  }
}

