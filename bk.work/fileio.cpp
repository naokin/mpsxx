#include <sstream>
#include "fileio.h"

std::string mpsxx::getfile (const std::string& name, const std::string& prefix, int i, int j, int k, int l)
{
  std::stringstream filename;
  filename << prefix << "/" << name << "." << i << "." << j << "." << k << "." << l << ".tmp";
  return filename.str();
}
