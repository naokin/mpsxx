#ifndef _MPSXX_CXX11_FERMIONIC_QUANTUM_H
#define _MPSXX_CXX11_FERMIONIC_QUANTUM_H 1

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

namespace mpsxx     {

namespace fermionic {

class Quantum {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_nptcl & m_nspin; }

public:
  const static Quantum zero() { return Quantum(0, 0); }

  Quantum() : m_nptcl(0), m_nspin(0) { }
  Quantum(int nptcl, int nspin) : m_nptcl(nptcl), m_nspin(nspin) { }

  bool operator== (const Quantum& other) const { return (m_nptcl == other.m_nptcl && m_nspin == other.m_nspin); }
  bool operator!= (const Quantum& other) const { return (m_nptcl != other.m_nptcl || m_nspin != other.m_nspin); }
  bool operator<  (const Quantum& other) const { return (m_nptcl == other.m_nptcl) ? (m_nspin < other.m_nspin) : (m_nptcl < other.m_nptcl); }
  bool operator>  (const Quantum& other) const { return (m_nptcl == other.m_nptcl) ? (m_nspin > other.m_nspin) : (m_nptcl > other.m_nptcl); }

  Quantum operator* (const Quantum& other) const { return Quantum(m_nptcl+other.m_nptcl, m_nspin+other.m_nspin); }

  Quantum operator+ () const { return Quantum(+m_nptcl, +m_nspin); }
  Quantum operator- () const { return Quantum(-m_nptcl, -m_nspin); }

  bool   parity () const { return m_nptcl & 1; }
  double clebsch() const { return 1.0; }

  const int& p() const { return m_nptcl; }
  const int& s() const { return m_nspin; }

  friend std::ostream& operator<< (std::ostream& ost, const Quantum& q) { return ost << "{" << std::setw(2) << q.m_nptcl << "," << std::setw(2) << q.m_nspin << "}"; }

private:
  int
    m_nptcl;
  int
    m_nspin;
};

}; // fermionic

}; // mpsxx

#endif // _MPSXX_CXX11_FERMIONIC_QUANTUM_H
