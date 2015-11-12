#ifndef _MPSXX_CXX11_PARTICLE_SYMMETRY_H
#define _MPSXX_CXX11_PARTICLE_SYMMETRY_H 1

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

namespace mpsxx     {

namespace fermionic {

class Particle {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_nptcl; }

public:
  const static Particle zero() { return Particle(0); }

  Particle() : m_nptcl(0) { }
  Particle(int nptcl) : m_nptcl(nptcl) { }

  bool operator== (const Particle& other) const { return m_nptcl == other.m_nptcl; }
  bool operator!= (const Particle& other) const { return m_nptcl != other.m_nptcl; }
  bool operator<  (const Particle& other) const { return m_nptcl <  other.m_nptcl; }
  bool operator>  (const Particle& other) const { return m_nptcl >  other.m_nptcl; }

  Particle operator* (const Particle& other) const { return Particle(m_nptcl+other.m_nptcl); }

  Particle operator+ () const { return Particle(+m_nptcl); }
  Particle operator- () const { return Particle(-m_nptcl); }

  bool   parity () const { return m_nptcl & 1; }
  double clebsch() const { return 1.0; }

  const int& p() const { return m_nptcl; }

  std::ostream& operator<< (std::ostream& ost, const Particle& q) { return ost << std::setw(2) << q.m_nptcl; }

private:
  int
    m_nptcl;
};

}; // fermionic

namespace bosonic {

class Particle {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_nptcl; }

public:
  const static Particle zero() { return Particle(0); }

  Particle() : m_nptcl(0) { }
  Particle(int nptcl) : m_nptcl(nptcl) { }

  bool operator== (const Particle& other) const { return m_nptcl == other.m_nptcl; }
  bool operator!= (const Particle& other) const { return m_nptcl != other.m_nptcl; }
  bool operator<  (const Particle& other) const { return m_nptcl <  other.m_nptcl; }
  bool operator>  (const Particle& other) const { return m_nptcl >  other.m_nptcl; }

  Particle operator* (const Particle& other) const { return Particle(m_nptcl+other.m_nptcl); }

  Particle operator+ () const { return Particle(+m_nptcl); }
  Particle operator- () const { return Particle(-m_nptcl); }

  bool   parity () const { return   0; }
  double clebsch() const { return 1.0; }

  const int& p() const { return m_nptcl; }

  std::ostream& operator<< (std::ostream& ost, const Particle& q) { return ost << std::setw(2) << q.m_nptcl; }

private:
  int
    m_nptcl;
};

}; // bosonic

}; // mpsxx

#endif // _MPSXX_CXX11_PARTICLE_SYMMETRY_H
