#ifndef __BTAS_SPARSE_VECTOR_HPP
#define __BTAS_SPARSE_VECTOR_HPP 1

#include <vector>
#include <pair>
#include <algorithm>

namespace btas {

/// pair comparison
template<typename _T1, typename _T2>
struct __less_comp {
   typedef std::pair<_T1, _T2> __arg_type;
   bool operator() (const __arg_type& x, const __arg_type& y) const { return x.first < y.first; }
};

/// sparse vector class
template<typename _T, class _Compare = __less_comp<unsigned long, _T>>
class SpVector {

   //
   //  Type names ==================================================================================
   //

   /// type of value
   typedef _T value_type;

   /// type of size
   typedef unsigned long size_type;

   /// type of compare function
   typedef _Compare compare;

   /// type of container
   typedef std::vector<std::pair<size_type, value_type>> container;

   /// iterator
   typedef typename container::iterator iterator;

   /// iterator to const
   typedef typename container::const_iterator const_iterator;

   //
   //  Constructors ================================================================================
   //

   /// default constructor
   SpVector()
   { }

   /// destructor
   ~SpVector()
   { }

   /// construct by sizes
   /// \param nmax total size of vector (in dense view)
   /// \param nnz number of estimated non-zero elements
   explicit
   SpVector(size_type nmax, size_type nnz = 0)
   : max_size_ (nmax)
   {
      if(nnz > 0) data_.reserve(nnz);
   }

   /// copy constructor
   SpVector(const SpVector& x)
   : max_size_ (x.max_size_)
   {
      size_type nnz = x.data_.size();
      data_.resize(nnz);
      NumericType<value_type>::fast_copy(nnz, x.data_.data(), 1, data_.data(), 1);
   }

   /// move constructor implemented in terms of move constructors of member
   SpVector(SpVector&& x)
   : max_size_ (x.max_size_), data_ (x.data_)
   { }

   //
   //  Assignment operators ========================================================================
   //

   /// copy assignment operator
   SpVector& operator= (const SpVector& x)
   {
      size_type nnz = x.data_.size();
      resize(x.max_size_, nnz);
      data_.resize(nnz);
      NumericType<value_type>::fast_copy(nnz, x.data_.data(), 1, data_.data(), 1);
      return *this;
   }

   /// move assignment operator implemented in terms of move assignments of member
   SpVector& operator= (SpVector&& x)
   {
      max_size_ = x.max_size_;
      data_ = x.data_;
      return *this;
   }

   //
   //  Size functions ==============================================================================
   //

   /// \return total size of vector
   size_type size() const { return max_size_; }

   /// \return number of non-zero elements
   size_type nnz() const { return data_.size(); }

   /// resize
   void resize(size_type nmax, size_type nnz = 0)
   {
      max_size_ = nmax;
      data_.clear();
      if(nnz > data_.capacity()) {
         data_.reserve(nnz);
      }
      else if(nnz > 0) {
         container().swap(data_); // reallocation
         data_.reserve(nnz);
      }
      // if(nnz == 0) capacity size is kept
   }

   //
   //  Iterators ===================================================================================
   //

   /// \return iterator to the first
   iterator begin() { return data_.begin(); }

   /// \return const iterator to the first
   const_iterator begin() const { return data_.begin(); }

   /// \return iterator to the end
   iterator end() { return data_.end(); }

   /// \return const iterator to the end
   const_iterator end() const { return data_.end(); }

   /// \return iterator to the n-th element
   iterator find(size_type n) { return std::lower_bound(data_.begin(), data_.end(), n, compare()); }

   /// \return const iterator to the n-th element
   const_iterator find(size_type n) const { return std::lower_bound(data_.begin(), data_.end(), n, compare()); }

   //
   //  Element access ==============================================================================
   //

   /// element access with range check
   /// if element doesn't exist, it's to be allocated
   value_type& operator[] (size_type n)
   {
      iterator it = this->find(n);
      if(it == data_.end())
      {
         MPSXX_RUNTIME_ASSERT(n < max_size_, "position is out of range");
         it = data_.insert(it, std::make_pair(n, NumericType<value_type>::zero()));
      }
      return it->second;
   }

   /// element access with range check
   const value_type& operator[] (size_type n) const
   {
      const_iterator it = this->find(n);
      MPSXX_RUNTIME_ASSERT(it != data_.end(), "required element doesn't exist");
      return it->second;
   }

   //
   //  Others ======================================================================================
   //

   /// sort elements by compare
   void sort()
   {
      std::sort(data_.begin(), data_.end(), compare());
   }

   /// test empty
   bool empty() const { return data_.empty(); }

   /// swap
   void swap(SpVector& x)
   {
      std::swap(max_size_, x.max_size_);
      data_.swap(x.data_); // reallocates memory
   }

   /// clear
   void clear()
   {
      max_size_ = 0;
      data_.clear(); // memory space is not deleted
   }

private:

   //
   //  Member variables ============================================================================
   //

   /// total size of vector
   size_type max_size_;

   /// non-zero elements stored as 1D-array
   container data_;

};

} // namespace btas

#endif // __BTAS_SPARSE_VECTOR_HPP
