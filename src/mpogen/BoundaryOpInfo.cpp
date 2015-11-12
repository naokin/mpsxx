#include <cassert>
#include "BoundaryOpInfo.h"

void mpsxx::BoundaryOpInfo::reset (size_t index)
{
  comp_ = DOT;
  data_.clear();

  size_t Nop = 0;

  // I
  data_.insert(std::make_pair(DMRG_OpType(op::QC::I),Nop++));
  // H
  data_.insert(std::make_pair(DMRG_OpType(op::QC::H,index),Nop++));

  // Cre
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB,index),Nop++));

  // Des
  data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB,index),Nop++));

  // CreComp <= CreCreDes
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_COMP,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_COMP,index),Nop++));

  // DesComp <= CreDesDes
  data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_COMP,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_COMP,index),Nop++));

  // CreCre, CreDes, DesCre, and DesDes
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_DESA,index,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_DESB,index,index),Nop++));

  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_CREB,index,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_DESB,index,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_CREB,index,index),Nop++));
  data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_DESB,index,index),Nop++));

  // check duplication
  assert(data_.size() == Nop);

  return;
}

void mpsxx::BoundaryOpInfo::reset (size_t L, size_t N, bool enable_swapswp)
{
  size_t R = N-L;
  size_t Nestm; // estimate # operators

  std::vector<size_t> loop_indices;

  if(!enable_swapswp || L < R) { /* loop indices = 0...L-1 */
    comp_ = RIGHT;
    for(size_t i = 0; i < L; ++i) loop_indices.push_back(i);
    // Estimated size of operators
    Nestm = 2+4*N+2*L+4*L*(2*L-1);
  }
  else      { /* loop indices = L...N-1 */
    comp_ = LEFT;
    for(size_t i = L; i < N; ++i) loop_indices.push_back(i);
    // Estimated size of operators
    Nestm = 2+4*N+2*R+4*R*(2*R-1);
  }
  data_.clear();

  size_t Nop = 0;

  // I
  if(enable_swapswp || R > 0)
    data_.insert(std::make_pair(DMRG_OpType(op::QC::I),Nop++));
  // H
  if((!enable_swapswp || R > 0) && L > 0)
    data_.insert(std::make_pair(DMRG_OpType(op::QC::H),Nop++));

  // case where the first & the last boundaries
  if(L == 0 || R == 0) return;

  if(enable_swapswp && L >= R) {
    for(size_t i = 0; i < L; ++i) {
      // Cre
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_COMP,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_COMP,i),Nop++));
      // Des
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_COMP,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_COMP,i),Nop++));
    }
    for(size_t i = L; i < N; ++i) {
      // CreComp
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB,i),Nop++));
      // DesComp
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB,i),Nop++));
    }
  }
  else {
    for(size_t i = 0; i < L; ++i) {
      // Cre
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB,i),Nop++));
      // Des
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB,i),Nop++));
    }
    for(size_t i = L; i < N; ++i) {
      // CreComp
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_COMP,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_COMP,i),Nop++));
      // DesComp
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_COMP,i),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_COMP,i),Nop++));
    }
  }

  // CreCre, CreDes, DesCre, and DesDes
  for(size_t i = 0; i < loop_indices.size(); ++i) {

    size_t ix = loop_indices[i];

    data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_DESA,ix,ix),Nop++));
    data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_DESB,ix,ix),Nop++));

    data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_CREB,ix,ix),Nop++));
    data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_DESB,ix,ix),Nop++));
    data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_CREB,ix,ix),Nop++));
    data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_DESB,ix,ix),Nop++));

    for(size_t j = 0; j < i; ++j) {

      size_t jx = loop_indices[j];

      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_CREA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_CREB,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_CREA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_CREB,ix,jx),Nop++));

      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_DESA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREA_DESB,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_DESA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::CREB_DESB,ix,jx),Nop++));

      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_CREA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_CREB,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_CREA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_CREB,ix,jx),Nop++));

      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_DESA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESA_DESB,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_DESA,ix,jx),Nop++));
      data_.insert(std::make_pair(DMRG_OpType(op::QC::DESB_DESB,ix,jx),Nop++));
    }
  }

  // check duplication
  assert(data_.size() == Nop);

  return;
}

std::vector<fermion> mpsxx::BoundaryOpInfo::qshape () const
{
  size_t Nop = data_.size();

  std::vector<fermion> qs(Nop);
  for(const_iterator it = data_.begin(); it != data_.end(); ++it)
    qs.at(it->second) = it->first.q();

  return qs;
}

void mpsxx::BoundaryOpInfo::clean (const std::vector<size_t>& ds)
{
  std::map<size_t,DMRG_OpType> invmap;

  for(auto it = data_.begin(); it != data_.end(); ++it) {
    if(ds[it->second] > 0)
      invmap.insert(std::make_pair(it->second,it->first));
  }

  data_.clear();

  size_t nnz = 0;
  for(auto it = invmap.begin(); it != invmap.end(); ++it, ++nnz)
    data_.insert(std::make_pair(it->second,nnz));
}

std::ostream& operator<< (std::ostream& ost, const mpsxx::BoundaryOpInfo& info)
{
  for(auto it = info.begin(); it != info.end(); ++it)
    ost << "\t[ " << it->second << " ] : " << it->first.label() << std::endl;
  return ost;
}
