
/*! \file  remove_index.h
 *  \brief external function to remove zero-sized indices in btas::QSDArray
 */
#ifndef _MPSXX_CXX11_REMOVE_INDEX_H
#define _MPSXX_CXX11_REMOVE_INDEX_H 1

#include <btas/QSPARSE/QSTArray.h>

namespace mpsxx {

template<size_t N, class Q>
void remove_index
(btas::QSTArray<double,N,Q>& x, const size_t& x_rank, btas::Dshapes& y_dn_shape)
{
  btas::Dshapes x_dn_shape = x.dshape(x_rank);
  // checking
  assert(x_dn_shape.size() == y_dn_shape.size());
  std::map<size_t, size_t> _index_map;
  size_t nnz = 0;
  for(size_t i = 0; i < x_dn_shape.size(); ++i) {
    if(x_dn_shape[i] > 0 && y_dn_shape[i] > 0) {
      assert(x_dn_shape[i] == y_dn_shape[i]);
      _index_map.insert(std::make_pair(i, nnz++));
    }
    else {
      y_dn_shape[i] = 0;
    }
  }
  btas::TVector<btas::Qshapes<Q>, N> x_q_shape = x.qshape();
  for(auto it = _index_map.begin(); it != _index_map.end(); ++it) {
    x_qi.push_back(x_q_shape[x_rank][it->first]);
  }
  x_q_shape[x_rank] = x_qi;
  btas::QSDArray x_new(x.q(), x_q_shape);
  for(auto ix = x.begin(); ix != x.end(); ++ix) {
    btas::IVector<N> _b_index = x.index(ix->first);
    auto it = _index_map.find(_b_index[x_rank]);
    if(it == _index_map.end()) continue;
    _b_index[rank] = it->second;
    x_new.insert(_b_index, *(ix->second));
  }
  x = std::move(x_new);
}

template<size_t M, size_t N, class Q>
void clean_index
(btas::QSTArray<double,M,Q>& x, const size_t& x_rank,
 btas::QSTArray<double,N,Q>& y, const size_t& y_rank)
{
}

};

#endif // _MPSXX_CXX11_REMOVE_INDEX_H
