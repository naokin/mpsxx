CXX=/home/naokin/gnu/gcc/4.7.1/bin/g++ -std=c++0x -I. -I/home/naokin/gnu/gcc/4.7.1/include -I/home/naokin/gnu/gcc/4.7.1/c++/4.7.1
CXXFLAGS=-g -O3 -fopenmp -D_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

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

SRC_DMRG = dmrg.C main.C btas_template_specialize.C ../driver/fileio.C

OBJ_DMRG = $(SRC_DMRG:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: dmrg

dmrg	: $(OBJ_DMRG)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o dmrg.x $(OBJ_DMRG) $(BTASLIB)

clean	:
	rm *.o; rm *.x; rm *.a;
