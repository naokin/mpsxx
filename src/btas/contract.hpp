


template<typename _T, unsigned long L, unsigned long M, unsigned long N, class Q>
void contract(
   const _T& alpha,
   const TArray<_T, L, Q>& A, const std::string& labelA,
   const TArray<_T, M, Q>& B, const std::string& labelB,
   const _T& beta,
         TArray<_T, N, Q>& C, const std::string& labelC)
{
   contract(alpha, A, Index<L>(labelA), B, Index<M>(labelB), beta, C, Index<N>(labelC));
}

template<typename _T, unsigned long L, unsigned long M, unsigned long N, class Q>
void contract(
   const _T& alpha,
   const TArray<_T, L, Q>& A, const Index<L>& indexA,
   const TArray<_T, M, Q>& B, const Index<M>& indexB,
   const _T& beta,
         TArray<_T, N, Q>& C, const Index<N>& indexC)
{
}
