#ifndef __MPSXX_MPOGEN_BOUNDARY_OPINFO_H
#define __MPSXX_MPOGEN_BOUNDARY_OPINFO_H 1

#include <iostream>
#include <map>
#include <cstdlib>

#include <legacy/QSPARSE/Qshapes.h>

#include "fermion.h"
#include "bit_operator_type.h"

namespace mpsxx {

//! DMRG operators at certain boundary
class boundary_opinfo {
public:
  typedef std::map<mpogen::BIT_OPERATOR_TYPE, size_t>::iterator       iterator;
  typedef std::map<mpogen::BIT_OPERATOR_TYPE, size_t>::const_iterator const_iterator;

  //! Direction
  enum DIRECTION {
    BACKWARD = 0, //!< left block >  right block
    FORWARD  = 1, //!< left block <= right block
    NOBOUND  = 2  //!< for single site
  };

  //! Default constructor
  boundary_opinfo() : m_bn_dir(NOBOUND) { }
  //! Constructor with site index
  explicit boundary_opinfo(size_t index) { reset(index); }
  //! Constructor with size
  boundary_opinfo(size_t L, size_t N, bool _enable_swap_sweep_dir = true) { reset(L, N, _enable_swap_sweep_dir); }
  //! Copy constructor
  boundary_opinfo(const boundary_opinfo& x) : m_bn_ops(x.m_bn_ops), m_bn_dir(x.m_bn_dir) { }

  //! Reset operators for single site
  void reset(size_t index);
  //! Reset operators
  void reset(size_t L, size_t N, bool _enable_swap_sweep_dir = true);

  //! Get/Set direction
  inline const DIRECTION& direction() const { return m_bn_dir; }
  inline       DIRECTION& direction()       { return m_bn_dir; }

  //! Number of operators
  inline size_t size() const { return m_bn_ops.size(); }

  //  Iterator functions

  inline const_iterator begin() const { return m_bn_ops.begin(); }
  inline       iterator begin()       { return m_bn_ops.begin(); }
  inline const_iterator end  () const { return m_bn_ops.end  (); }
  inline       iterator end  ()       { return m_bn_ops.end  (); }

  inline const_iterator find(mpogen::BIT_OPERATOR_TYPE _type) const { return m_bn_ops.find(_type); }
  inline       iterator find(mpogen::BIT_OPERATOR_TYPE _type)       { return m_bn_ops.find(_type); }

  btas::Qshapes<fermion> get_qshape() const;

  void clean(const btas::Dshapes& _dn_shape);

private:
  DIRECTION
    m_bn_dir;
  std::map<mpogen::BIT_OPERATOR_TYPE, size_t>
    m_bn_ops;
};

} // mpsxx

std::ostream& operator<< (std::ostream& ost, const mpsxx::boundary_opinfo& info);

#endif // __MPSXX_MPOGEN_BOUNDARY_OPINFO_H
