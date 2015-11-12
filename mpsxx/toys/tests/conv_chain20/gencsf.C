#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
using namespace std;

int main()
{
  int L = 20;

  vector<int> index(L, 0);

  int iend = 0;
  int icsf = 0;
  while(!iend) {
    int nspin = accumulate(index.begin(), index.end(), 0);
    if(nspin == L/2) {
      cout << "\tCSF[" << setw(8) << icsf++ << "]: ";
      for(int i = 0; i < L; ++i) {
        if(index[i] == 0) cout << "d";
        else              cout << "u";
      }
      cout << endl;
    }
    for(int i = L-1; i >= 0; --i) {
      if(++index[i] < 2) break;
      index[i] = 0;
      if(i == 0) iend = 1;
    }
  }

  return 0;
}
