all: testfourier

GSL_INC = /lab/software/apparatus3/cpptools/gsl-1.15/
APP3_CPP_INC = /lab/software/apparatus3/cpp/
CPP_TOOLS_INC = /lab/software/apparatus3/cpptools/
INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC}

GSL_LIB = -L/lab/software/apparatus3/cpptools/gsl-1.15/.libs/ -L/lab/software/apparatus3/cpptools/gsl-1.15/cblas/.libs
CCFITS_LIB = /lab/software/apparatus3/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /lab/software/apparatus3/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/lab/software/apparatus3/cpptools/gsl-1.15/.libs/:/lab/software/apparatus3/cpptools/gsl-1.15/cblas/.libs/:/lab/software/apparatus3/cpptools/CCfits/.libs/:/lab/software/apparatus3/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS = -Wall ${INC} 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm  ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS}

objs =  /lab/software/apparatus3/cpp/utils/utils.o  /lab/software/apparatus3/cpp/qini/qini_utils.o

testfourier: testfourier.o fourier.o
	g++  $^ ${objs} ${LFLAGS} -o testfourier

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




