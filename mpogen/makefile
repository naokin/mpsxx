CXX=icpc -std=c++0x
CXXFLAGS=-g -O3 -openmp -D_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=
BLASINC=
BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lguide -lpthread

BOOSTDIR=/home/naokin/boost/1.54.0
BOOSTINC=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASDIR=/home/naokin/btas-0.3.0/c++11
BTASINC=-I. -I.. -I$(BTASDIR)/include
BTASLIB= $(BTASDIR)/lib/libbtas.a

INCLUDEFLAGS= $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

SRC_MPOGEN  = bit_operator_type.C boundary_opinfo.C generate_qc_operators.C generate_site_operator.C get_product_ops.C integral_component.C parsing_fcidump.C prime_operators.C mpogen.C ../driver/fileio.C
SRC_HUBGEN  = bit_operator_type.C boundary_opinfo.C generate_hubbard_operators.C generate_site_operator.C get_product_ops.C integral_component.C parsing_fcidump.C prime_operators.C hubgen.C ../driver/fileio.C
SRC_TESTGEN = bit_operator_type.C boundary_opinfo.C get_product_ops.C testgen.C

OBJ_MPOGEN  = $(SRC_MPOGEN:.C=.o)
OBJ_HUBGEN  = $(SRC_HUBGEN:.C=.o)
OBJ_TESTGEN = $(SRC_TESTGEN:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: mpogen hubgen testgen

mpogen	: $(OBJ_MPOGEN)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o mpogen.x $(OBJ_MPOGEN) $(BTASLIB)

hubgen	: $(OBJ_HUBGEN)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o hubgen.x $(OBJ_HUBGEN) $(BTASLIB)

testgen	: $(OBJ_TESTGEN)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o testgen.x $(OBJ_TESTGEN) $(BTASLIB)

clean	:
	rm *.o; rm *.x; rm *.a;
