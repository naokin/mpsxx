#ifndef __MPSXX_MATRIX_STATE_HPP
#define __MPSXX_MATRIX_STATE_HPP 1

namespace mpsxx {

template<typename _T, class _Quantum>
class MatrixState {
public:
   typedef _T value_type;
   typedef _Quantum qnum_type;
   typedef btas::QSTArray<value_type, 2, qnum_type> matrix_type;
   typedef boost::shared_ptr<matrix_type> matrix_ptr;
private:
   std::vector<boost::shared_ptr<matrix_type>> store_;
};

} // namespace mpsxx

#endif // __MPSXX_MATRIX_STATE_HPP
