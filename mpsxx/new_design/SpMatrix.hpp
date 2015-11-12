#ifndef __MPSXX_SPARSE_MATRIX_HPP
#define __MPSXX_SPARSE_MATRIX_HPP 1

namespace mpsxx {

/// iterator
template<class _Key, class _T>
class __matrix_iterator : public typename std::map<_Key, _T>::iterator {

   typedef typename std::map<_Key, _T>::iterator iterator;

public:

   __matrix_iterator() : stride_(1) { }

   __matrix_iterator(iterator pos, size_type n) : iterator(pos), stride_(n) { }

   /// copy constructor
   __matrix_iterator(const __matrix_iterator& i) : iterator(i), stride_(i.stride_) { }

   /// copy assignment
   __matrix_iterator& operator= (const __matrix_iterator& i) {
      iterator::operator=(i);
      stride_ = i.stride_;
      return *this;
   }

   size_type row() const { return this->first/stride_; }

   size_type col() const { return this->first%stride_; }

private:

   size_type stride_;

};

/// iterator
template<class _Key, class _T>
class __matrix_const_iterator : public typename std::map<_Key, _T>::const_iterator {

   typedef typename std::map<_Key, _T>::const_iterator iterator;

public:

   __matrix_const_iterator() : stride_(1) { }

   __matrix_const_iterator(iterator pos, size_type n) : iterator(pos), stride_(n) { }

   /// copy constructor
   __matrix_const_iterator(const __matrix_const_iterator& i) : iterator(i), stride_(i.stride_) { }

   /// copy assignment
   __matrix_const_iterator& operator= (const __matrix_const_iterator& i) {
      iterator::operator=(i);
      stride_ = i.stride_;
      return *this;
   }

   size_type row() const { return this->first/stride_; }

   size_type col() const { return this->first%stride_; }

private:

   size_type stride_;

};

template<typename _T>
class SpMatrix {

public:

   /// element type
   typedef _T value_type;

   /// size type
   typedef unsigned long size_type;

   /// storage type
   typedef std::map<size_type, value_type> storage_type;

   /// iterator
   typedef __matrix_iterator<size_type, value_type> iterator;

   /// const iterator
   typedef __matrix_const_iterator<size_type, value_type> const_iterator;

   /// default constructor
   SpMatrix() : nrows_(0), ncols_(0) { }

   /// construct by sizes
   SpMatrix(size_type nr, size_type nc) : nrows_(nr), ncols_(nc) { }

   /// copy constructor
   SpMatrix(const SpMatrix& x) : nrows_(x.nrows_), ncols_(x.ncols_), data_(x.data_) { }

   /// element access w/ allocation
   value_type& operator() (size_type i, size_type j) {
      assert(i < nrows_ && j < ncols_);
      // map::insert returns pair<iterator, bool> in C++11
      // if key has already existed, this inserts nothing to return iterator to existed value
      return *(data_.insert(std::make_pair(i*nrows_+j, value_type())).first);
   }

   /// element access
   const value_type& operator() (size_type i, size_type j) const {
      assert(i < nrows_ && j < ncols_);
      auto it = data_.find(i*nrows_+j);
      assert(it != data_.end());
      return *it;
   }

   //
   // size information
   //

   size_type size() const { return nrows_*ncols_; }

   size_type nnz() const { return data_.size(); }

   const size_type& rows() const { return nrows_; }

   const size_type& cols() const { return ncols_; }

   //
   // get index
   //

   size_type  rows(const iterator& it) const { return it->first/ncols_; }

   size_type  rows(const const_iterator& it) const { return it->first/ncols_; }

   size_type  cols(const iterator& it) const { return it->first%ncols_; }

   size_type  cols(const const_iterator& it) const { return it->first%ncols_; }

   //
   // return iterator
   //

   iterator begin() { return data_.begin(); }

   const_iterator begin() const { return data_.begin(); }

   iterator end() { return data_.end(); }

   const_iterator end() const { return data_.end(); }

private:

   /// # rows
   size_type nrows_;

   /// # cols
   size_type ncols_;

   /// map to elements
   storage_type data_;

};

};

#endif // __MPSXX_SPARSE_MATRIX_HPP
