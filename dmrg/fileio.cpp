#include <sstream>

#include "mpidefs.h"

#include "fileio.h"

std::string mpsxx::getfile (const std::string& name, const std::string& prefix, int site, int iState, int jState)
{
  Communicator world;
  std::stringstream filename;
  filename << prefix << "/" << name << "." << site << "." << iState << "." << jState << "." << world.rank() << ".tmp";
  return filename.str();
}
