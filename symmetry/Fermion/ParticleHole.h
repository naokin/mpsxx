#ifndef _MPSXX_CXX11_PARTICLE_HOLE_SYMMETRY_H
#define _MPSXX_CXX11_PARTICLE_HOLE_SYMMETRY_H 1

namespace mpsxx     {

namespace fermionic {

class ParticleHole {
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) { ar & m_np_ref & m_ph_sym; }

public:
  const static ParticleHole zero() { return ParticleHole(0, 0); }

  ParticleHole(int _np_ref = 0, int _ph_lev = 0) : m_np_ref(_np_ref), m_ph_lev(_ph_lev) { }

  bool operator== (const ParticleHole& other) const { return (m_np_ref == other.m_np_ref && m_ph_lev == other.m_ph_lev); }
  bool operator!= (const ParticleHole& other) const { return (m_np_ref != other.m_np_ref || m_ph_lev != other.m_ph_lev); }
  bool operator<  (const ParticleHole& other) const { return (m_np_ref == other.m_np_ref) ? (m_ph_lev < other.m_ph_lev) : (m_np_ref < other.m_np_ref); }
  bool operator>  (const ParticleHole& other) const { return (m_np_ref == other.m_np_ref) ? (m_ph_lev > other.m_ph_lev) : (m_np_ref > other.m_np_ref); }

  ParticleHole operator* (const ParticleHole& other) const { return ParticleHole(m_np_ref+other.m_np_ref, m_ph_lev+other.m_ph_lev); }

  ParticleHole operator+ (const ParticleHole& other) const { return ParticleHole(m_np_ref+other.m_np_ref, m_ph_lev+other.m_ph_lev); }
  ParticleHole operator- (const ParticleHole& other) const { return ParticleHole(m_np_ref-other.m_np_ref, m_ph_lev-other.m_ph_lev); }

  ParticleHole operator+ () const { return ParticleHole(+m_np_ref, +m_ph_lev); }
  ParticleHole operator- () const { return ParticleHole(-m_np_ref, -m_ph_lev); }

  const int& p() const { return m_np_ref+m_ph_lev; }

  const int& np_ref() const { return m_np_ref; }
  const int& ph_lev() const { return m_ph_lev; }

private:
  //! Number of particles in ref. state
  int m_np_ref;

  //! Number of particles(+)/holes(-) added to ref. state
  int m_ph_lev;
};

}; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_PARTICLE_HOLE_SYMMETRY_H
