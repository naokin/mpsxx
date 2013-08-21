#include <iostream>

#include "Fermion/ParticleHole.h"
#include "Spin/Spin.h"
#include "PointGroup/D2h.h"
#include "GenericQuantum.h"

using namespace mpsxx;

int main() {
  typedef fermionic::ParticleHole PQ;
  typedef Spin                    SQ;
  typedef pointgroup::D2h         PG;

  typedef GenericQuantum<PQ, SQ, PG> MyQuantum;

  MyQuantum q1(PQ(1,1,0), SQ(+1), PG(PG::B1u));
  MyQuantum q2(PQ(1,0,1), SQ(-1), PG(PG::B3g));
  MyQuantum q3 = MyQuantum::zero();

  std::cout << "q1 = " << q1 << std::endl;
//std::cout << "q2 = " << q2 << std::endl;
//std::cout << "q3 = " << q3 << std::endl;

//if(q1 == q2) std::cout << "q1 == q2" << std::endl;
//if(q1 != q2) std::cout << "q1 != q2" << std::endl;
//if(q1 <  q2) std::cout << "q1 <  q2" << std::endl;
//if(q1 >  q2) std::cout << "q1 >  q2" << std::endl;

//q3 = q1 * q2;
//std::cout << "q3 := q1 * q2 = " << q3 << std::endl;

  return 0;
}
