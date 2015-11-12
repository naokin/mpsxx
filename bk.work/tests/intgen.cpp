#include <iostream>
#include <iomanip>
#include <fstream>

int main ()
{
  using namespace std;

  size_t N; cin >> N;

  ofstream fo("FCIDUMP");

  fo << " &FCI NORB=" << setw(3) << N << ",NELEC=" << setw(3) << N << ",MS2= 0," << endl;
  fo << "  ORBSYM="; for(size_t i = 0; i < N; ++i) fo << "1,"; fo << endl;
  fo << "  ISYM=1," << endl;
  fo << " /" << endl;

  bool FTii   = false;
  bool FTij   = true;
  bool FViiii = true;
  bool FViiij = false;
  bool FViijj = false;
  bool FVijij = false;
  bool FVijjj = false;
  bool FViijk = false;
  bool FVijik = false;
  bool FVijjk = false;
  bool FVjjik = false;
  bool FVijkk = false;
  bool FVikjk = false;
  bool FVijkl = true;
  bool FVikjl = true;
  bool FViljk = true;

  fo.precision(16);

  for(size_t i = 0; i < N; ++i) {
    if(FViiii) {
      fo << setw(24) << scientific << 1.0 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << i+1 << endl;
    }
    for(size_t j = 0; j < i; ++j) {
      if(FViiij) {
        fo << setw(24) << scientific << 0.5 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << j+1 << endl;
      }
      if(FViijj) {
        fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << j+1 << endl;
      }
      if(FVijij) {
        fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << i+1 << setw(4) << j+1 << endl;
      }
      if(FVijjj) {
        fo << setw(24) << scientific << 0.5 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << j+1 << setw(4) << j+1 << endl;
      }
      for(size_t k = 0; k < j; ++k) {
        if(FViijk) {
          fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
        }
        if(FVijik) {
          fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << i+1 << setw(4) << k+1 << endl;
        }
        if(FVijjk) {
          fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
        }
        if(FVjjik) {
          fo << setw(24) << scientific << 0.2 << setw(4) << j+1 << setw(4) << j+1 << setw(4) << i+1 << setw(4) << k+1 << endl;
        }
        if(FVijkk) {
          fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << setw(4) << k+1 << endl;
        }
        if(FVikjk) {
          fo << setw(24) << scientific << 0.2 << setw(4) << i+1 << setw(4) << k+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
        }
        for(size_t l = 0; l < k; ++l) {
          if(FVijkl) {
            fo << setw(24) << scientific << 0.1 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << setw(4) << l+1 << endl;
          }
          if(FVikjl) {
            fo << setw(24) << scientific << 0.1 << setw(4) << i+1 << setw(4) << k+1 << setw(4) << j+1 << setw(4) << l+1 << endl;
          }
          if(FViljk) {
            fo << setw(24) << scientific << 0.1 << setw(4) << i+1 << setw(4) << l+1 << setw(4) << j+1 << setw(4) << k+1 << endl;
          }
        }
      }
    }
  }
  for(size_t i = 0; i < N; ++i) {
    if(FTii) {
      fo << setw(24) << scientific << -1.0 << setw(4) << i+1 << setw(4) << i+1 << setw(4) << 0 << setw(4) << 0 << endl;
    }
    if(FTij) {
      for(size_t j = 0; j < i; ++j) {
        fo << setw(24) << scientific << 1.0 << setw(4) << i+1 << setw(4) << j+1 << setw(4) << 0 << setw(4) << 0 << endl;
      }
    }
  }
  fo << setw(24) << scientific << 0.0 << setw(4) << 0 << setw(4) << 0 << setw(4) << 0 << setw(4) << 0 << endl;

  return 0;
}
