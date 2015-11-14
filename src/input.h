#ifndef __MPSXX_DMRG_INPUT_H
#define __MPSXX_DMRG_INPUT_H

#include <vector>
#include <string>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

namespace mpsxx {

enum DMRGALGORITHM { ONESITE, TWOSITE };

struct DMRGInput {

  bool
    restart;
  std::vector<double>
    energy;

  size_t
    N_sites;

  int
    N_spins;
  int
    N_elecs;

  size_t
    N_roots;
  int
    N_max_states;

  DMRGALGORITHM
    algorithm;

  size_t
    N_max_sweep_iter;
  double
    tolerance;
  double
    noise;

  std::string
    prefix;

  DMRGInput()
  : restart          ( 0 ),
    N_sites          ( 0 ),
    N_spins          ( 0 ),
    N_elecs          ( 0 ),
    N_roots          ( 1 ),
    N_max_states     ( 0 ),
    algorithm        ( ONESITE ),
    N_max_sweep_iter ( 100 ),
    tolerance        ( 1.0e-8 ),
    noise            ( 0.0 ),
    prefix           ( "." )
  { }

  DMRGInput (const std::string& fname)
  : DMRGInput () // Deligation (C++11)
  { this->parse(fname); }

  void parse (const std::string& fname);

private:

  friend class boost::serialization::access;

  /// Boost serialization
  template<class Archive>
  void serialize (Archive& ar, const unsigned int version)
  {
    ar & restart;
    ar & energy;
    ar & N_sites;
    ar & N_spins;
    ar & N_elecs;
    ar & N_roots;
    ar & N_max_states;
    ar & algorithm;
    ar & N_max_sweep_iter;
    ar & tolerance;
    ar & noise;
    ar & prefix;
  }

};

} // namespace mpsxx

#endif // __MPSXX_DMRG_INPUT_H
