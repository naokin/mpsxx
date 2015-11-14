Tensor<double,N,CblasRowMajor>;

SpTensor<double,N,CblasRowMajor>;

QSpTensor<double,N,Q,CblasRowMajor>;

QSpMerge<N,Q,CblasRowMajor>; // ?

{
  QSpTensor<double,5,Q,CblasRowMajor> t(...);
  QSpTensor<double,3,Q,CblasRowMajor> m;

  QSpMerge<5,3,Q,CblasRowMajor> info(t,shape(0,1),shape(2),shape(3,4));

  m = t.merge(info);
  t = m.expand(info);
}
