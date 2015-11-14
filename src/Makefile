CXX=mpic++
CXXFLAGS=-g -std=c++11 -O3 -D_DEFAULT_QUANTUM -D_ENABLE_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=/home100/opt/intel/mkl
BLASINC=-I$(BLASDIR)/include
BLASLIB=-L$(BLASDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

BOOSTINC=-I/home100/naokin/boost/1.59.0-gcc/include
BOOSTLIB=-L/home100/naokin/boost/1.59.0-gcc/lib -lboost_serialization -lboost_filesystem -lboost_mpi

BTASINC=-I../include

INCLUDEFLAGS= -I. -I.. $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

SRC_MPOGEN  = bit_operator_type.cpp boundary_opinfo.cpp gen_qc_operators.cpp gen_site_operator.cpp get_product_ops.cpp integral_component.cpp parsing_integral.cpp prime_operators.cpp make_qc_mpos.cpp mpogen.cpp fileio.cpp compress_qc_mpos.cpp
SRC_NAIVE_MPOGEN  = gen_qc_naive_mpos.cpp compress_qc_mpos.cpp parsing_integral.cpp naive_mpogen.cpp fileio.cpp
SRC_DMRG = dmrg.cpp make_sweep.cpp make_random_mpss.cpp fileio.cpp input.cpp
SRC_UTIL = opex.cpp fileio.cpp

OBJ_MPOGEN  = $(SRC_MPOGEN:.cpp=.o)
OBJ_NAIVE_MPOGEN  = $(SRC_NAIVE_MPOGEN:.cpp=.o)
OBJ_DMRG = $(SRC_DMRG:.cpp=.o)
OBJ_UTIL = $(SRC_UTIL:.cpp=.o)

.cpp.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: mpogen naive_mpogen dmrg opex

mpogen	: $(OBJ_MPOGEN)
	$(CXX) $(CXXFLAGS) -o mpogen.x $(OBJ_MPOGEN) $(LIBRARYFLAGS)

naive_mpogen	: $(OBJ_NAIVE_MPOGEN)
	$(CXX) $(CXXFLAGS) -o naive_mpogen.x $(OBJ_NAIVE_MPOGEN) $(LIBRARYFLAGS)

dmrg	: $(OBJ_DMRG)
	$(CXX) $(CXXFLAGS) -o dmrg.x $(OBJ_DMRG) $(LIBRARYFLAGS)

opex	: $(OBJ_UTIL)
	$(CXX) $(CXXFLAGS) -o opex.x $(OBJ_UTIL) $(LIBRARYFLAGS)

clean	:
	rm *.o; rm *.x; rm *.a;
