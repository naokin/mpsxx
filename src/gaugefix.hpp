#ifndef __MPSXX_GAUGE_FIX_HPP
#define __MPSXX_GAUGE_FIX_HPP

#include <iostream>
#include <iomanip>

#include <btas/QSPARSE/QSTArray.h>

#include "input.h"
#include "canonicalize.hpp"
#include "compute_guesswave.hpp"
#include "fileio.h"

namespace mpsxx {

/// Compute expectation value <i|O|j>
template<class Q>
void gaugefix (const DMRGInput& input, const size_t& iroot = 0)
{
  using std::cout;
  using std::endl;
  using std::setw;

  size_t N = input.N_sites;
  size_t K = input.N_roots;

  assert(iroot < K);

  btas::QSTArray<double,3,Q> lmps;
  btas::QSTArray<double,3,Q> rmps;
  btas::QSTArray<double,3,Q> wave;

  if(0) { // algo 1

  load(wave,getfile("wave",input.prefix,0,iroot));
  for(size_t i = 0; i < N-1; ++i) {
    load(rmps,getfile("rmps",input.prefix,i+1,iroot));
    canonicalize(1,wave,lmps,rmps);
    save(lmps,getfile("lmps",input.prefix,i,iroot));
    wave = rmps;
  }
  for(size_t i = N-1; i > 0; --i) {
    btas::Normalize(wave);
    save(wave,getfile("wave",input.prefix,i,iroot));

    load(lmps,getfile("lmps",input.prefix,i-1,iroot));
    canonicalize(0,wave,lmps,rmps);
    save(rmps,getfile("rmps",input.prefix,i,iroot));
    wave = lmps;
  }
    btas::Normalize(wave);
    save(wave,getfile("wave",input.prefix,0,iroot));

  } else { // algo 2

  load(wave,getfile("wave",input.prefix,0,iroot));
  for(size_t i = 0; i < N-1; ++i) {
    load(lmps,getfile("lmps",input.prefix,i,  iroot));
    load(rmps,getfile("rmps",input.prefix,i+1,iroot));
    btas::QSTArray<double,3,Q> temp;
    compute_guesswave(1,lmps,wave,rmps,temp);
    wave = temp;
    save(wave,getfile("wave",input.prefix,i+1,iroot));
  }
  for(size_t i = N-1; i > 0; --i) {
    load(rmps,getfile("rmps",input.prefix,i,  iroot));
    load(lmps,getfile("lmps",input.prefix,i-1,iroot));
    btas::QSTArray<double,3,Q> temp;
    compute_guesswave(0,rmps,wave,lmps,temp);
    wave = temp;
    save(wave,getfile("wave",input.prefix,i-1,iroot));
  }

  } // end if
}

} // namespace mpsxx

#endif // __MPSXX_GAUGE_FIX_HPP
