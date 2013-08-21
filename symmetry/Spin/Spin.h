
#ifndef _MPSXX_CXX11_SPIN_QUANTUM_H
#define _MPSXX_CXX11_SPIN_QUANTUM_H 1

#include <iostream>
#include <boost/serialization/serialization.hpp>

namespace mpsxx {

//! Spin quantum number class
class Spin {
private:
  //! Boost serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_spin; }

public:
  const static Spin zero() { return Spin(0); }

  Spin() : m_spin(0) { }
  Spin(int spin) : m_spin(spin) { }

  inline bool operator== (const Spin& other) const { return m_spin == other.m_spin; }
  inline bool operator!= (const Spin& other) const { return m_spin != other.m_spin; }
  inline bool operator<  (const Spin& other) const { return m_spin <  other.m_spin; }
  inline bool operator>  (const Spin& other) const { return m_spin >  other.m_spin; }

  inline Spin operator* (const Spin& other) const { return Spin(m_spin+other.m_spin); }
  inline Spin operator+ () const { return Spin(+m_spin); }
  inline Spin operator- () const { return Spin(-m_spin); }

  const int& s() const { return m_spin; }

  friend std::ostream& operator<< (std::ostream& ost, const Spin& q) { return ost << q.m_spin; }

private:
  //! Spin quantum number
  int m_spin;
};

}; // namespace mpsxx

#endif // _MPSXX_CXX11_SPIN_QUANTUM_H
