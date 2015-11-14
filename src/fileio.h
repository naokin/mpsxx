#ifndef __MPSXX_FILEIO_H
#define __MPSXX_FILEIO_H

#include <fstream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

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

/// Get the filename "[prefix]/[name].i.j.k.l"
std::string getfile(const std::string& name, const std::string& prefix, int site = 0, int iState = 0, int jState = 0);

} // namespace mpsxx

#endif // __MPSXX_FILEIO_H
