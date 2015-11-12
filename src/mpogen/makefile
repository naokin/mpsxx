CXX=g++ -std=c++11
CXXFLAGS=-g -O2 -D_DEFAULT_QUANTUM -D_ENABLE_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=/home100/opt/intel/mkl
BLASINC=-I$(BLASDIR)/include
BLASLIB=-L$(BLASDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

BOOSTINC=-I/homec/naokin/boost/1.54.0/include
BOOSTLIB=-L/homec/naokin/boost/1.54.0/lib -lboost_serialization

BTASDIR= $(HOME)/devel/btas
BTASINC= -I$(BTASDIR)/include
BTASLIB= $(BTASDIR)/lib/libbtas.a

INCLUDEFLAGS= -I. -I.. $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

#SRC_MPOGEN  = add_operators.cpp BoundaryOpInfo.cpp Category.cpp gen_product_operators.cpp gen_qcmpo.cpp gen_site_operators.cpp Get1eComponent.cpp Get2eComponent.cpp
SRC_MPOGEN  = add_operators.cpp BoundaryOpInfo.cpp Category.cpp gen_product_operators.cpp Get1eComponent.cpp Get2eComponent.cpp test.cpp

OBJ_MPOGEN  = $(SRC_MPOGEN:.cpp=.o) ../driver/fileio.o ../driver/parsing_integral.o

.cpp.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: mpogen

mpogen	: $(OBJ_MPOGEN)
	$(CXX) $(CXXFLAGS) -o mpogen.x $(OBJ_MPOGEN) $(BTASLIB) $(LIBRARYFLAGS)

test	: $(OBJ_MPOGEN)
	$(CXX) $(CXXFLAGS) -o test.x $(OBJ_MPOGEN) $(BTASLIB) $(LIBRARYFLAGS)

clean	:
	rm *.o; rm *.x; rm *.a;
