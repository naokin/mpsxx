#ifndef _MPSXX_CXX11_FILEIO_H
#define _MPSXX_CXX11_FILEIO_H 1

#include <fstream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <MPX_Types.h>

template<class T>
void load(      T& obj, const std::string& f_name)
{
  std::ifstream f_load(f_name.c_str());
  boost::archive::binary_iarchive ia(f_load);
  ia & obj;
  return;
}

template<class T>
void save(const T& obj, const std::string& f_name)
{
  std::ofstream f_save(f_name.c_str());
  boost::archive::binary_oarchive oa(f_save);
  oa & obj;
  return;
}

namespace mpsxx {

enum MPS_TYPE { WAVEFUNCTION, LEFTCANONICAL, RIGHTCANONICAL };

//std::string get_mpofile(const std::string& prefix, const MPO_TYPE& _type, const int& index);
std::string get_mpofile(const std::string& prefix,                        const int& index);
std::string get_mpsfile(const std::string& prefix, const MPS_TYPE& _type, const int& index);
std::string get_oprfile(const std::string& prefix, const MPS_TYPE& _type, const int& index);

}; // namespace mpsxx

#endif // _MPSXX_CXX11_FILEIO_H
