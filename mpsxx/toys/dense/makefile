CXX=g++
CXXFLAGS=-g -O3 -fopenmp

BLASDIR=/home100/opt/intel/mkl
BLASINC=-I$(BLASDIR)/include
BLASLIB=-L$(BLASDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

BOOSTDIR=/homec/naokin/boost/1.54.0
BOOSTINC=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASDIR=/homec/naokin/devel/btas
BTASINC=-I$(BTASDIR)/include
BTASLIB=  $(BTASDIR)/lib

INCLUDEFLAGS=-I. $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS=    $(BLASLIB) $(BOOSTLIB)

SRC_SAMPLE = main.C dmrg.C lrt.C driver.C btas_template_specialize.C

OBJ_SAMPLE = $(SRC_SAMPLE:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: dmrg.x

dmrg.x	: $(OBJ_SAMPLE)
	$(CXX) $(CXXFLAGS) $(LIBRARYFLAGS) -o dmrg.x $(OBJ_SAMPLE) $(BTASLIB)/libbtas.a

clean	:
	rm *.o; rm *.x; rm *.a;
