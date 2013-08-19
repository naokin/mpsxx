#ifndef _MPSXX_FERMIONIC_QUANTUM_H
#define _MPSXX_FERMIONIC_QUANTUM_H 1

#include <iostream>
#include <iomanip>
#include <boost/serialization/serialization.hpp>

namespace mpsxx     {

namespace fermionic {

class Quantum
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_np_ref & m_ph_sym; }

public:
  const static Quantum zero() { return Quantum(0, 0, 0); }

  Quantum(int _np_ref = 0, int _ph_lev = 0, int _n_spin = 0) : m_np_ref(_np_ref), m_ph_lev(_ph_lev), m_n_spin(_n_spin) { }

  bool operator== (const Quantum& other) const { return (m_np_ref == other.m_np_ref && m_ph_lev == other.m_ph_lev && m_n_spin == other.m_n_spin); }

  bool operator!= (const Quantum& other) const { return (m_np_ref != other.m_np_ref || m_ph_lev != other.m_ph_lev || m_n_spin != other.m_n_spin); }

  bool operator<  (const Quantum& other) const {
    if(m_np_ref != other.m_np_ref) return m_np_ref < other.m_np_ref;
    if(m_ph_lev != other.m_ph_lev) return m_ph_lev < other.m_ph_lev;
                                   return m_n_spin < other.m_n_spin;
  }

  bool operator>  (const Quantum& other) const {
    if(m_np_ref != other.m_np_ref) return m_np_ref > other.m_np_ref;
    if(m_ph_lev != other.m_ph_lev) return m_ph_lev > other.m_ph_lev;
                                   return m_n_spin > other.m_n_spin;
  }

  Quantum operator* (const Quantum& other) const { return std::move(Quantum(m_np_ref+other.m_np_ref, m_ph_lev+other.m_ph_lev, m_n_spin+other.m_n_spin)); }

  Quantum operator+ (const Quantum& other) const { return std::move(Quantum(m_np_ref+other.m_np_ref, m_ph_lev+other.m_ph_lev, m_n_spin+other.m_n_spin)); }
  Quantum operator- (const Quantum& other) const { return std::move(Quantum(m_np_ref+other.m_np_ref, m_ph_lev+other.m_ph_lev, m_n_spin+other.m_n_spin)); }

  Quantum operator+ () const { return std::move(Quantum(+m_np_ref, +m_ph_lev, +m_n_spin)); }
  Quantum operator- () const { return std::move(Quantum(-m_np_ref, -m_ph_lev, -m_n_spin)); }

  bool parity() const { return (m_np_ref+m_ph_lev) & 1; }

  const int& p() const { return m_np_ref+m_ph_lev; }
  const int& s() const { return m_n_spin; }

  const int& np_ref() const { return m_np_ref; }
  const int& ph_lev() const { return m_ph_lev; }

private:
  //! Number of particles in ref. state
  int m_np_ref;

  //! Number of particles(+)/holes(-) added to ref. state
  int m_ph_lev;

  //! Number of spins
  int m_n_spin;
};

}; // namespace fermionic

}; // namespace mpsxx

inline std::ostream& operator<< (std::ostream& ost, const mpsxx::fermionic::Quantum& q) {
  ost << "(" << std::setw(2) << q.p() << "[" << std::setw(2) << m_ph_lev << "]," << std::setw(2) << q.s() << ")";
  return ost;
}

#endif // _MPSXX_FERMIONIC_QUANTUM_H
