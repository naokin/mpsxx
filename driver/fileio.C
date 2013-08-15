
#include <driver/fileio.h>

//std::string mpsxx::get_mpofile
//(const std::string& prefix, const mpsxx::MPO_TYPE& _type, const int& index)
//{
//  std::stringstream filename;
//  switch(_type) {
//    case HEISENBERG:
//      filename << prefix << "heisenberg-mpo_site-" << index << /* mpigetrank() << */ ".tmp";
//      break;
//    case FERMION_HUBBARD:
//      filename << prefix << "fermion-hubbard-mpo_site-" << index << /* mpigetrank() << */ ".tmp";
//      break;
//    case MOLECULAR:
//      filename << prefix << "qc-mpo_site-" << index << /* mpigetrank() << */ ".tmp";
//      break;
//    default:
//      break;
//  }
//  return filename.str();
//}

std::string mpsxx::get_mpofile
(const std::string& prefix, const int& index)
{
  std::stringstream filename;
  filename << prefix << "mpo_site-" << index << /* mpigetrank() << */ ".tmp";
  return filename.str();
}

std::string mpsxx::get_mpsfile
(const std::string& prefix, const mpsxx::MPS_TYPE& _type, const int& index)
{
  std::stringstream filename;
  switch(_type) {
    case WAVEFUNCTION:
      filename << prefix << "wave-mps_site-" << index << /* mpigetrank() << */ ".tmp";
      break;
    case LEFTCANONICAL:
      filename << prefix << "left-mps_site-" << index << /* mpigetrank() << */ ".tmp";
      break;
    case RIGHTCANONICAL:
      filename << prefix << "right-mps_site-" << index << /* mpigetrank() << */ ".tmp";
      break;
    default:
      break;
  }
  return filename.str();
}

std::string mpsxx::get_oprfile
(const std::string& prefix, const mpsxx::MPS_TYPE& _type, const int& index)
{
  std::stringstream filename;
  switch(_type) {
//  case WAVEFUNCTION:
//    filename << prefix << "trans-opr_site-" << index << /* mpigetrank() << */ ".tmp";
//    break;
    case LEFTCANONICAL:
      filename << prefix << "left-opr_site-" << index << /* mpigetrank() << */ ".tmp";
      break;
    case RIGHTCANONICAL:
      filename << prefix << "right-opr_site-" << index << /* mpigetrank() << */ ".tmp";
      break;
    default:
      break;
  }
  return filename.str();
}

