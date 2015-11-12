#ifndef __MPSXX_MPOGEN_BOUNDARY_OP_INFO_H
#define __MPSXX_MPOGEN_BOUNDARY_OP_INFO_H

#include <ostream>
#include <vector>
#include <map>

#include <symmetry/fermion.h>
#include "DMRG_OpType.h"

namespace mpsxx {

/// DMRG normal & complementary operators for a certain boundary
struct BoundaryOpInfo {

  typedef std::map<DMRG_OpType,size_t> data_t;

  typedef data_t::iterator iterator;

  typedef data_t::const_iterator const_iterator;

  /// of which block is a complementary block?
  enum COMPLEMENTARY_BLOCK {
    DOT   = 0, ///< dot block
    LEFT  = 1, ///< in case sizeof(left) >  sizeof(right)
    RIGHT = 2  ///< in case sizeof(left) <= sizeof(right)
  };

  /// Default constructor
  BoundaryOpInfo () : comp_(DOT) { }

  /// Constructor with site index
  explicit BoundaryOpInfo (size_t index) { this->reset(index); }

  /// Constructor with size
  BoundaryOpInfo (size_t L, size_t N, bool enable_swapswp = true) { this->reset(L,N,enable_swapswp); }

  /// Copy constructor
  BoundaryOpInfo(const BoundaryOpInfo& x) : comp_(x.comp_), data_(x.data_) { }

  /// Reset operators for single site
  void reset (size_t index);

  /// Reset operators
  void reset (size_t L, size_t N, bool enable_swapswp = true);

  /// Get flag of complementary block
  inline const COMPLEMENTARY_BLOCK& comp () const { return comp_; }

  /// Number of operators
  inline size_t size () const { return data_.size(); }

  //  Iterator functions

  inline const_iterator begin() const { return data_.begin(); }
  inline       iterator begin()       { return data_.begin(); }
  inline const_iterator end  () const { return data_.end  (); }
  inline       iterator end  ()       { return data_.end  (); }

  inline const_iterator find(const DMRG_OpType& type) const { return data_.find(type); }
  inline       iterator find(const DMRG_OpType& type)       { return data_.find(type); }

  std::vector<fermion> qshape() const;

  /// remove zero block by given size array \c ds for this boundary
  void clean (const std::vector<size_t>& ds);

private:

  COMPLEMENTARY_BLOCK comp_;

  data_t data_;

};

} // mpsxx

std::ostream& operator<< (std::ostream& ost, const mpsxx::BoundaryOpInfo& info);

#endif // __MPSXX_MPOGEN_BOUNDARY_OP_INFO_H
