
MatrixState::
TArray<SpMatrix<Matrix<double>>,1> // rank 3

MatrixOperator::
TArray<SpMatrix<Matrix<double>>,2> // rank 4

BlockOperator::
SpVector<Matrix<double>>

SVD::
TArray<SpMatrix<Matrix<double>>,2> // rank 4

/// forward SVD for one-site MPS
SVD (const TArray<SpMatrix<Matrix<double>>,1>& A,
     SpVector<Vector<double>>& S, TArray<SpMatrix<Matrix<double>>,1>& U, SpMatrix<Matrix<double>>& V)

/// backward SVD for one-site MPS
SVD (const TArray<SpMatrix<Matrix<double>>,1>& A,
     SpVector<Vector<double>>& S, SpMatrix<Matrix<double>>& U, TArray<SpMatrix<Matrix<double>>,1>& V)

/// SVD for two-site MPS
SVD (const QTArray<QSpMatrix<Matrix<double>>,2>& A,
     SpVector<Vector<double>>& S, QTArray<QSpMatrix<Matrix<double>>,1>& U, QTArray<QSpMatrix<Matrix<double>>,1>& V)
{
   Matrix<double> Ax;
}

MpArray<double, 3, 1, 0, 2> // MPS

MpArray<double, 4, 2, 2, 0> // MPO

MpArray<double, 3, 0, 1, 2> // BOP

/// Array container class for MPS language
/// \tparam T type of value
/// \tparam N rank of object on top
/// \tparam Phys rank of physical index
/// \tparam Sp rank of sparse object
/// \tparam Dn rank of dense object
template<typename T, unsigned int N, class... IndexGroups>
class MpArray {
};
