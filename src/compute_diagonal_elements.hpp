#ifndef __MPSXX_COMPUTE_DIAGONAL_ELEMENTS_H
#define __MPSXX_COMPUTE_DIAGONAL_ELEMENTS_H

#include <btas/QSPARSE/QSTArray.h>

#include "mpidefs.h"

namespace mpsxx {

//! Compute diagonal H elements for one-site algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSTArray<double,4,Q>& mpo0,
 const btas::QSTArray<double,3,Q>& lopr,
 const btas::QSTArray<double,3,Q>& ropr,
       btas::QSTArray<double,3,Q>& diag)
{
  btas::STArray<double,3> mpo0_diag;
  btas::STArray<double,2> lopr_diag;
  btas::STArray<double,2> ropr_diag;

  btas::Tie(mpo0,btas::shape(1,2),mpo0_diag);
  btas::Tie(lopr,btas::shape(0,2),lopr_diag);
  btas::Tie(ropr,btas::shape(0,2),ropr_diag);

  btas::STArray<double,3> scr1;
  btas::Contract(1.0,lopr_diag,btas::shape(1),mpo0_diag,btas::shape(0),1.0,scr1);
  btas::STArray<double,3> scr2;
  btas::Contract(1.0,scr1,     btas::shape(2),ropr_diag,btas::shape(1),1.0,scr2);

  btas::Copy(scr2,diag,1); // preserve quantum number of diag
}

//! Compute diagonal H elements for two-site algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSTArray<double,4,Q>& lmpo,
 const btas::QSTArray<double,4,Q>& rmpo,
 const btas::QSTArray<double,3,Q>& lopr,
 const btas::QSTArray<double,3,Q>& ropr,
       btas::QSTArray<double,4,Q>& diag)
{
  btas::STArray<double,3> lmpo_diag;
  btas::STArray<double,3> rmpo_diag;
  btas::STArray<double,2> lopr_diag;
  btas::STArray<double,2> ropr_diag;

  btas::Tie(lmpo,btas::shape(1,2),lmpo_diag);
  btas::Tie(rmpo,btas::shape(1,2),rmpo_diag);
  btas::Tie(lopr,btas::shape(0,2),lopr_diag);
  btas::Tie(ropr,btas::shape(0,2),ropr_diag);

  btas::STArray<double,3> scr1;
  btas::Contract(1.0,lopr_diag,btas::shape(1),lmpo_diag,btas::shape(0),1.0,scr1);
  btas::STArray<double,4> scr2;
  btas::Contract(1.0,scr1,     btas::shape(2),rmpo_diag,btas::shape(0),1.0,scr2);
  btas::STArray<double,4> scr3;
  btas::Contract(1.0,scr2,     btas::shape(3),ropr_diag,btas::shape(1),1.0,scr3);

  btas::Copy(scr3,diag,1); // preserve quantum number of diag
}

//! Compute diagonal H elements for merged algorithm
template<class Q>
void compute_diagonal_elements
(const btas::QSTArray<double,3,Q>& lopr,
 const btas::QSTArray<double,3,Q>& ropr,
       btas::QSTArray<double,2,Q>& diag)
{
  btas::STArray<double,2> lopr_diag;
  btas::STArray<double,2> ropr_diag;

  btas::Tie(lopr,btas::shape(0,2),lopr_diag);
  btas::Tie(ropr,btas::shape(0,2),ropr_diag);

  btas::STArray<double,2> scr1;
  btas::Gemm(btas::NoTrans,btas::Trans,1.0,lopr_diag,ropr_diag,1.0,scr1);

  btas::Copy(scr1,diag,1); // preserve quantum number of diag
}

//! Compute diagonal H elements for one-site algorithm
template<class Q>
void compute_diagonal_elements
(const std::vector<btas::QSTArray<double,4,Q>>& mpo0,
 const std::vector<btas::QSTArray<double,3,Q>>& lopr,
 const std::vector<btas::QSTArray<double,3,Q>>& ropr,
       btas::QSTArray<double,3,Q>& diag)
{
  assert(mpo0.size() == lopr.size());
  assert(mpo0.size() == ropr.size());
  for(size_t g = 0; g < mpo0.size(); ++g) {
    btas::QSTArray<double,3,Q> diag_tmp(diag.q(),diag.qshape());
    compute_diagonal_elements(mpo0[g],lopr[g],ropr[g],diag_tmp);
    btas::Axpy(1.0,diag_tmp,diag);
    diag.check_dshape();
  }
#ifndef _SERIAL
  Communicator world;
  for(size_t i = 1; i < world.size(); ++i) {
    if(world.rank() == i) {
      world.send(0,i,diag);
    }
    if(world.rank() == 0) {
      btas::QSTArray<double,3,Q> temp;
      world.recv(i,i,temp);
      btas::Axpy(1.0,temp,diag);
      diag.check_dshape();
    }
  }
  boost::mpi::broadcast(world,diag,0);
#endif
}

//! Compute diagonal H elements for two-site algorithm
template<class Q>
void compute_diagonal_elements
(const std::vector<btas::QSTArray<double,4,Q>>& lmpo,
 const std::vector<btas::QSTArray<double,4,Q>>& rmpo,
 const std::vector<btas::QSTArray<double,3,Q>>& lopr,
 const std::vector<btas::QSTArray<double,3,Q>>& ropr,
       btas::QSTArray<double,4,Q>& diag)
{
  assert(lmpo.size() == rmpo.size());
  assert(lmpo.size() == lopr.size());
  assert(rmpo.size() == ropr.size());
  for(size_t g = 0; g < lmpo.size(); ++g) {
    btas::QSTArray<double,4,Q> diag_tmp(diag.q(),diag.qshape());
    compute_diagonal_elements(lmpo[g],rmpo[g],lopr[g],ropr[g],diag_tmp);
    btas::Axpy(1.0,diag_tmp,diag);
    diag.check_dshape();
  }
#ifndef _SERIAL
  Communicator world;
  for(size_t i = 1; i < world.size(); ++i) {
    if(world.rank() == i) {
      world.send(0,i,diag);
    }
    if(world.rank() == 0) {
      btas::QSTArray<double,4,Q> temp;
      world.recv(i,i,temp);
      btas::Axpy(1.0,temp,diag);
      diag.check_dshape();
    }
  }
  boost::mpi::broadcast(world,diag,0);
#endif
}

//! Compute diagonal H elements for merged algorithm
template<class Q>
void compute_diagonal_elements
(const std::vector<btas::QSTArray<double,3,Q>>& lopr,
 const std::vector<btas::QSTArray<double,3,Q>>& ropr,
       btas::QSTArray<double,2,Q>& diag)
{
  assert(lopr.size() == ropr.size());
  for(size_t g = 0; g < lopr.size(); ++g) {
    btas::QSTArray<double,2,Q> diag_tmp(diag.q(),diag.qshape());
    compute_diagonal_elements(lopr[g],ropr[g],diag_tmp);
    btas::Axpy(1.0,diag_tmp,diag);
    diag.check_dshape();
  }
#ifndef _SERIAL
  Communicator world;
  for(size_t i = 1; i < world.size(); ++i) {
    if(world.rank() == i) {
      world.send(0,i,diag);
    }
    if(world.rank() == 0) {
      btas::QSTArray<double,2,Q> temp;
      world.recv(i,i,temp);
      btas::Axpy(1.0,temp,diag);
      diag.check_dshape();
    }
  }
  boost::mpi::broadcast(world,diag,0);
#endif
}

} // namespace mpsxx

#endif // __MPSXX_COMPUTE_DIAGONAL_ELEMENTS_H
