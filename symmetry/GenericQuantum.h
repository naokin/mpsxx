#ifndef _MPSXX_CXX11_GENERIC_QUANTUM_H
#define _MPSXX_CXX11_GENERIC_QUANTUM_H 1

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>

namespace mpsxx {

//! Generic quantum number container class
/*!
 *  \param Q  first quantum number class
 *  \param Qs second and later quantum number classes
 *
 *  Note that at least 1 quantum number class must be assigned
 *  Varidable template expression enables to describe a set of quantum number classes
 *
 *  e.g.)
 *
 *  GenericQuantum<fermionic::Particle, Spin, D2h>
 *
 *  gives a set of quantum numbers of electronic system with D2h point group symmetry
 */
template<class Q, class... Qs>
class GenericQuantum : public Q, public Qs... {
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    boost_serialize::mf_serialize<Archive, Q, Qs...>(ar, version);
  }

  template<class Archive>
  struct boost_serialize {
    template<class P, class... Ps>
    inline void mf_serialize(Archive ar, const unsigned int version) { ar & boost::base_object<P>(*this); mf_serialize<Archive, Ps...>(ar, version); }
    inline void mf_serialize(Archive ar, const unsigned int version) { }
  };

  template<class P, class... Ps>
  inline bool mf_equal_to    (const P& x, const Ps&... xs) const { return (static_cast<const P&>(*this) == x) ? mf_equal_to    (xs...) : false; }
  inline bool mf_equal_to    () const { return true ; }

  template<class P, class... Ps>
  inline bool mf_not_equal_to(const P& x, const Ps&... xs) const { return (static_cast<const P&>(*this) == x) ? mf_not_equal_to(xs...) : true ; }
  inline bool mf_not_equal_to() const { return false; }

  template<class P, class... Ps>
  inline bool mf_less        (const P& x, const Ps&... xs) const { return (static_cast<const P&>(*this) == x) ? mf_less        (xs...) : (static_cast<const P&>(*this) < x); }
  inline bool mf_less        () const { return false; }

  template<class P, class... Ps>
  inline bool mf_greater     (const P& x, const Ps&... xs) const { return (static_cast<const P&>(*this) == x) ? mf_greater     (xs...) : (static_cast<const P&>(*this) > x); }
  inline bool mf_greater     () const { return false; }

  template<class P, class... Ps>
  inline void mf_copy_assignment(const P & x, const Ps &... xs) { static_cast<P&>(*this) = x; mf_copy_assignment(xs...); }
  inline void mf_copy_assignment() { }

  template<class P, class... Ps>
  inline void mf_move_assignment(      P&& x,       Ps&&... xs) { static_cast<P&>(*this) = x; mf_move_assignment(xs...); }
  inline void mf_move_assignment() { }

//template<class P, class... Ps>
//inline double mf_clebsch(const P& x, const P&... xs) const { return static_cast<const P&>(x).clebsch() * mf_clebsch(xs...); }
//inline double mf_clebsch() const { return 1.0; }

  template<class P, class... Ps>
  inline std::ostream& mf_printf(std::ostream& ost, const P& x, const Ps&... xs) const { return ost << "," << x << mf_printf(ost, xs...); }
  inline std::ostream& mf_printf(std::ostream& ost) const { return ost; }

public:
  const static GenericQuantum zero() { return GenericQuantum(Q::zero(), Qs::zero()...); }

  GenericQuantum() : { }

  GenericQuantum(const Q& q, const Qs&... qs)  : Q(q), Qs(qs)... { }

  GenericQuantum(const GenericQuantum&  other) : Q(static_cast<const Q& >(other)), Qs(static_cast<const Qs& >(other))... { }
  GenericQuantum(      GenericQuantum&& other) : Q(static_cast<      Q&&>(other)), Qs(static_cast<      Qs&&>(other))... { }

  GenericQuantum& operator= (const GenericQuantum&  other) { mf_copy_assignment(static_cast<const Q &>(other), static_cast<const Qs &>(other)...); return *this; }
  GenericQuantum& operator= (      GenericQuantum&& other) { mf_copy_assignment(static_cast<      Q&&>(other), static_cast<      Qs&&>(other)...); return *this; }

  bool operator== (const GenericQuantum& other) const { return mf_equal_to    (static_cast<const Q&>(other), static_cast<const Qs&>(other)...); }
  bool operator!= (const GenericQuantum& other) const { return mf_not_equal_to(static_cast<const Q&>(other), static_cast<const Qs&>(other)...); }
  bool operator<  (const GenericQuantum& other) const { return mf_less        (static_cast<const Q&>(other), static_cast<const Qs&>(other)...); }
  bool operator>  (const GenericQuantum& other) const { return mf_greater     (static_cast<const Q&>(other), static_cast<const Qs&>(other)...); }

  GenericQuantum operator* (const GenericQuantum& other) const { return GenericQuantum(Q::operator* (static_cast<const Q&>(other)), Qs::operator* (static_cast<const Qs&>(other))...); }

  GenericQuantum operator+ () const { return GenericQuantum(Q::operator+ (), Qs::operator+ ()...); }
  GenericQuantum operator- () const { return GenericQuantum(Q::operator- (), Qs::operator- ()...); }

//double clebsch() const { return mf_clebsch(static_cast<const Q&>(*this), static_cast<const Qs&>(*this)...); }

  friend std::ostream& operator<< (std::ostream& ost, const GenericQuantum& q) { return ost << " { " << static_cast<const Q&>(q) << q.mf_printf(ost, static_cast<const Qs&>(q)...) << " } "; }

};

}; // namespace mpsxx

#endif // _MPSXX_CXX11_GENERIC_QUANTUM_H
