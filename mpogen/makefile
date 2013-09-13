CXX=g++ -std=c++11
CXXFLAGS=-g -O3 -openmp -D_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=
BLASINC=
BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential

BOOSTLIB=-lboost_serialization

BTASDIR=/home/bright/bestanden/programmas/btas
BTASINC=-I. -I.. -I$(BTASDIR)/include
BTASLIB= $(BTASDIR)/lib/libbtas.a

INCLUDEFLAGS= $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

SRC_MPOGEN  = bit_operator_type.C boundary_opinfo.C generate_qc_operators.C generate_site_operator.C get_product_ops.C integral_component.C ../driver/parsing_integral.C prime_operators.C mpogen.C ../driver/fileio.C
SRC_HUBGEN  = bit_operator_type.C boundary_opinfo.C generate_hubbard_operators.C generate_site_operator.C get_product_ops.C integral_component.C ../driver/parsing_integral.C prime_operators.C hubgen.C ../driver/fileio.C
SRC_TESTGEN = bit_operator_type.C boundary_opinfo.C get_product_ops.C testgen.C

OBJ_MPOGEN  = $(SRC_MPOGEN:.C=.o)
OBJ_HUBGEN  = $(SRC_HUBGEN:.C=.o)
OBJ_TESTGEN = $(SRC_TESTGEN:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: mpogen hubgen testgen

mpogen	: $(OBJ_MPOGEN)
	$(CXX) $(CXXFLAGS) -o mpogen.x $(OBJ_MPOGEN) $(BTASLIB) $(LIBRARYFLAGS)

hubgen	: $(OBJ_HUBGEN)
	$(CXX) $(CXXFLAGS) -o hubgen.x $(OBJ_HUBGEN) $(BTASLIB) $(LIBRARYFLAGS)

testgen	: $(OBJ_TESTGEN)
	$(CXX) $(CXXFLAGS) -o testgen.x $(OBJ_TESTGEN) $(BTASLIB) $(LIBRARYFLAGS)

clean	:
	rm *.o; rm *.x; rm *.a;
