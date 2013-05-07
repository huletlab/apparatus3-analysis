all: testfourier

GSL_INC = /home/ernie/analysis/cpptools/gsl-1.15/
APP3_CPP_INC = /home/ernie/analysis/apparatus3-analysis/
CPP_TOOLS_INC = /home/ernie/analysis/cpptools/
INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC}

GSL_LIB = -L/home/ernie/analysis/cpptools/gsl-1.15/.libs/ -L/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs
CCFITS_LIB = /home/ernie/analysis/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /home/ernie/analysis/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/home/ernie/analysis/cpptools/gsl-1.15/.libs/:/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs/:/home/ernie/analysis/cpptools/CCfits/.libs/:/home/ernie/analysis/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS = -Wall ${INC} 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm  ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS}

objs =  /home/ernie/analysis/apparatus3-analysis/utils/utils.o  /home/ernie/analysis/apparatus3-analysis/qini/qini_utils.o

testfourier: testfourier.o fourier.o
	g++  $^ ${objs} ${LFLAGS} -o testfourier

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




