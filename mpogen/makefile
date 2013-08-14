CXX=/home/naokin/gnu/gcc/4.7.1/bin/g++ -std=c++0x -I. -I/home/naokin/gnu/gcc/4.7.1/include -I/home/naokin/gnu/gcc/4.7.1/c++/4.7.1
CXXFLAGS=-g -O2 -fopenmp -D_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=
BLASINC=
BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lguide -lpthread

BOOSTDIR=/home/naokin/boost/1.52.0
BOOSTINC=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASDIR=/home/naokin/btas-0.3.0/c++11
BTASINC=-I. -I.. -I$(BTASDIR)/include
BTASLIB= $(BTASDIR)/lib/libbtas.a

INCLUDEFLAGS= $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

SRC_MPOGEN  = bit_operator_type.C boundary_opinfo.C generate_qc_operators.C generate_site_operator.C get_product_ops.C integral_component.C parsing_fcidump.C prime_operators.C mpogen.C ../driver/fileio.C
SRC_TESTGEN = bit_operator_type.C boundary_opinfo.C get_product_ops.C testgen.C

OBJ_MPOGEN  = $(SRC_MPOGEN:.C=.o)
OBJ_TESTGEN = $(SRC_TESTGEN:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: mpogen testgen

mpogen	: $(OBJ_MPOGEN)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o mpogen.x $(OBJ_MPOGEN) $(BTASLIB)

testgen	: $(OBJ_TESTGEN)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o testgen.x $(OBJ_TESTGEN) $(BTASLIB)

clean	:
	rm *.o; rm *.x; rm *.a;
