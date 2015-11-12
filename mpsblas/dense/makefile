CXX=/homec/naokin/gnu/gcc/4.8.4/bin/g++
CXXFLAGS=-g -std=c++11 -O3 -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=/home100/opt/intel/mkl
BLASINC=-I$(BLASDIR)/include
BLASLIB=-L$(BLASDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

BOOSTDIR=/homec/naokin/boost/1.54.0
BOOSTINC=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASDIR=/homec/naokin/devel/btas
BTASINC=-I$(BTASDIR)/include
BTASLIB=-L$(BTASDIR)/lib -lbtas

INCLUDEFLAGS=-I. $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS=    $(BLASLIB) $(BOOSTLIB) $(BTASLIB)

SRC_DMRG = main.cpp hubbard.cpp optimize.cpp
OBJ_DMRG = $(SRC_DMRG:.cpp=.o)

.cpp.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: dmrg.x

dmrg.x	: $(OBJ_DMRG)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o dmrg.x $(OBJ_DMRG)

clean	:
	rm *.o; rm *.x; rm *.a;
