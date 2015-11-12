#include <iostream>
#include <vector>

int main () {
  if(std::vector<double>(4).size()  == 4) std::cout << "OK" << std::endl;
  return 0;
}
