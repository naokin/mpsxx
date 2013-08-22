#ifndef _MPSXX_CXX11_PARTICLE_HOLE_SYMMETRY_H
#define _MPSXX_CXX11_PARTICLE_HOLE_SYMMETRY_H 1

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

namespace mpsxx     {

namespace fermionic {

class ParticleHole {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_np & m_pc & m_hc; }

public:
  const static ParticleHole zero() { return ParticleHole(0, 0, 0); }

  ParticleHole(int _np = 0, int _pc = 0, int _hc = 0) : m_np(_np), m_pc(_pc), m_hc(_hc) { }

  bool operator== (const ParticleHole& other) const { return (m_np == other.m_np && m_pc == other.m_pc && m_hc == other.m_hc); }
  bool operator!= (const ParticleHole& other) const { return (m_np != other.m_np || m_pc != other.m_pc || m_hc != other.m_hc); }

  bool operator<  (const ParticleHole& other) const {
    if(m_np == other.m_np) {
      if(m_pc == other.m_pc) {
        return (m_hc < other.m_hc);
      }
      else {
        return (m_pc < other.m_pc);
      }
    }
    else {
        return (m_np < other.m_np);
    }
  }
  bool operator>  (const ParticleHole& other) const {
    if(m_np == other.m_np) {
      if(m_pc == other.m_pc) {
        return (m_hc > other.m_hc);
      }
      else {
        return (m_pc > other.m_pc);
      }
    }
    else {
        return (m_np > other.m_np);
    }
  }

  ParticleHole operator* (const ParticleHole& other) const { return ParticleHole(m_np+other.m_np, m_pc+other.m_pc, m_hc+other.m_hc); }

  ParticleHole operator+ () const { return ParticleHole(+m_np, +m_pc, +m_hc); }
  ParticleHole operator- () const { return ParticleHole(-m_np, -m_pc, -m_hc); }

  bool   parity () const { return m_np & 1; }
  double clebsch() const { return 1.0; }

  const int& p () const { return m_np; }

  const int& pc() const { return m_pc; }
  const int& hc() const { return m_hc; }

  friend std::ostream& operator<< (std::ostream& ost, const ParticleHole& q) { return ost << std::setw(2) << q.m_np << ":[" << std::setw(2) << q.m_pc << "," << std::setw(2) << q.m_hc << "]"; }

private:
  //! Number of particles
  int m_np;

  //! Number of particles created to the ref. state
  int m_pc;

  //! Number of holes     created to the ref. state
  int m_hc;
};

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_PARTICLE_HOLE_SYMMETRY_H
