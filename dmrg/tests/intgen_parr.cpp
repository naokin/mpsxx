#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>

bool is_local (size_t nproc, size_t iproc, size_t i, size_t j)
{
  return (i % nproc == iproc);
}

bool is_local (size_t nproc, size_t iproc, size_t i, size_t j, size_t k, size_t l)
{
  size_t ij = i*(i+1)/2+j;
  return (ij % nproc == iproc);
}

int main ()
{
  using namespace std;

  size_t N; cin >> N;

  ofstream fo("FCIDUMP");

  fo << " &FCI NORB=" << setw(3) << N << ",NELEC=" << setw(3) << N << ",MS2= 0," << endl;
  fo << "  ORBSYM="; for(size_t i = 0; i < N; ++i) fo << "1,"; fo << endl;
  fo << "  ISYM=1," << endl;
  fo << " /" << endl;

  bool FTii   = true;
  bool FTij   = true;
  bool FViiii = true;
  bool FViiij = true;
  bool FViijj = true;
  bool FVijij = true;
  bool FVijjj = true;
  bool FViijk = true;
  bool FVijik = true;
  bool FVijjk = true;
  bool FVjjik = true;
  bool FVijkk = true;
  bool FVikjk = true;
  bool FVijkl = true;
  bool FVikjl = true;
  bool FViljk = true;

  mt19937 rgen; uniform_real_distribution<double> dist(0.0,1.0);

  fo.precision(16);

  size_t iproc = 0;
  size_t nproc = 1;

  for(size_t i = 0; i < N; ++i) {
    if(FViiii && is_local(nproc,iproc,i,i,i,i)) {
      fo << setw(24) << scientific << 1.0*dist(rgen) << setw(4) << i+1 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << i+1 << endl;
    }
    for(size_t j = 0; j < i; ++j) {
      if(FViiij && is_local(nproc,iproc,i,i,i,j)) {
        fo << setw(24) << scientific << 0.5*dist(rgen) << setw(4) << i+1 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << j+1 << endl;
      }
      if(FViijj && is_local(nproc,iproc,i,i,j,j)) {
        fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << j+1 << endl;
      }
      if(FVijij && is_local(nproc,iproc,i,i,j,j)) {
        fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << i+1 << setw(4) << j+1 << endl;
      }
      if(FVijjj && is_local(nproc,iproc,i,j,j,j)) {
        fo << setw(24) << scientific << 0.5*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << j+1 << setw(4) << j+1 << endl;
      }
      for(size_t k = 0; k < j; ++k) {
        if(FViijk && is_local(nproc,iproc,i,i,j,k)) {
          fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
        }
        if(FVijik && is_local(nproc,iproc,i,i,j,k)) {
          fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << i+1 << setw(4) << k+1 << endl;
        }
        if(FVijjk && is_local(nproc,iproc,i,j,j,k)) {
          fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
        }
        if(FVjjik && is_local(nproc,iproc,i,j,j,k)) {
          fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << j+1 << setw(4) << j+1 << setw(4) << i+1 << setw(4) << k+1 << endl;
        }
        if(FVijkk && is_local(nproc,iproc,i,j,k,k)) {
          fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << setw(4) << k+1 << endl;
        }
        if(FVikjk && is_local(nproc,iproc,i,j,k,k)) {
          fo << setw(24) << scientific << 0.2*dist(rgen) << setw(4) << i+1 << setw(4) << k+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
        }
        for(size_t l = 0; l < k; ++l) {
          if(FVijkl && is_local(nproc,iproc,i,j,k,l)) {
            fo << setw(24) << scientific << 0.1*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << setw(4) << l+1 << endl;
          }
          if(FVikjl && is_local(nproc,iproc,i,j,k,l)) {
            fo << setw(24) << scientific << 0.1*dist(rgen) << setw(4) << i+1 << setw(4) << k+1 << setw(4) << j+1 << setw(4) << l+1 << endl;
          }
          if(FViljk && is_local(nproc,iproc,i,j,k,l)) {
            fo << setw(24) << scientific << 0.1*dist(rgen) << setw(4) << i+1 << setw(4) << l+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
          }
        }
      }
    }
  }
  for(size_t i = 0; i < N; ++i) {
    if(FTii && is_local(nproc,iproc,i,i)) {
      fo << setw(24) << scientific << -1.0*dist(rgen) << setw(4) << i+1 << setw(4) << i+1 << setw(4) << 0 << setw(4) << 0 << endl;
    }
    for(size_t j = 0; j < i; ++j) {
      if(FTij && is_local(nproc,iproc,i,j)) {
        fo << setw(24) << scientific << 1.0*dist(rgen) << setw(4) << i+1 << setw(4) << j+1 << setw(4) << 0 << setw(4) << 0 << endl;
      }
    }
  }
  fo << setw(24) << scientific << 0.0 << setw(4) << 0 << setw(4) << 0 << setw(4) << 0 << setw(4) << 0 << endl;

  return 0;
}
