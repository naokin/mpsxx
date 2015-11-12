#ifndef __MPSXX_MPOGEN_DMRG_OPTYPE_H
#define __MPSXX_MPOGEN_DMRG_OPTYPE_H

#include <ostream>
#include <array>

#include "Category.h"
#include <symmetry/fermion.h>

namespace mpsxx {

struct DMRG_OpType
{
  op::QC::CATEGORY category;

  std::array<int,2> index;

  DMRG_OpType (op::QC::CATEGORY c, int i = -1, int j = -1)
  {
    category = c;
    index[0] = i;
    index[1] = j;
  }

  inline fermion q () const { return op::QC::QNumTable[category]; }

  inline std::string label () const { return op::QC::LabTable[category]; }

inline bool operator== (const DMRG_OpType& x) const
{
   if(this->category != x.category) return false;
   return (this->index[0] == x.index[0] && this->index[1] == x.index[1]);
}

inline bool operator!= (const DMRG_OpType& x) const
{
   if(this->category != x.category) return true;
   return (this->index[0] != x.index[0] || this->index[1] != x.index[1]);
}

inline bool operator< (const DMRG_OpType& x) const
{
   if(this->category != x.category) return this->category < x.category;
   return (this->index[0] != x.index[0]) ? (this->index[0] < x.index[0]) : (this->index[1] < x.index[1]);
}

inline bool operator<= (const DMRG_OpType& x) const
{
   if(this->category != x.category) return this->category < x.category;
   return (this->index[0] != x.index[0]) ? (this->index[0] < x.index[0]) : (this->index[1] <= x.index[1]);
}

inline bool operator> (const DMRG_OpType& x) const
{
   if(this->category != x.category) return this->category > x.category;
   return (this->index[0] != x.index[0]) ? (this->index[0] > x.index[0]) : (this->index[1] > x.index[1]);
}

inline bool operator>= (const DMRG_OpType& x) const
{
   if(this->category != x.category) return this->category > x.category;
   return (this->index[0] != x.index[0]) ? (this->index[0] > x.index[0]) : (this->index[1] >= x.index[1]);
}

};

} // namespace mpsxx

inline std::ostream& operator<< (std::ostream& ost, const mpsxx::DMRG_OpType& op)
{
  ost << op.label();
  if(op.index[0] >= 0) {
    ost << "(" << op.index[0];
    if(op.index[1] >= 0) {
      ost << "," << op.index[1];
    }
    ost << ")";
  }
  return ost;
}

#endif // __MPSXX_MPOGEN_DMRG_OPTYPE_H
