#ifndef _MPSXX_CXX11_H
#define _MPSXX_CXX11_H 1

#include <vector>
#include <cassert>

#include <iostream>
#define MPSXX_DEBUG(msg)\
{ std::cout << "MPSXX_DEBUG: " << msg << std::endl; }

#include <stdexcept>
#define MPSXX_THROW(truth, msg)\
{ if (!(truth)) { throw std::runtime_error(msg); } }

namespace mpsxx {

typedef unsigned int  uint;
typedef unsigned long wint;

}; // namespace btas

#endif // _MPSXX_CXX11_H
