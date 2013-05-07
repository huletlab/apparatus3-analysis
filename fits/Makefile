all: testfits

GSL_INC = /home/ernie/analysis/cpptools/gsl-1.15/
APP3_CPP_INC = /home/ernie/analysis/apparatus3-analysis/
CPP_TOOLS_INC = -I/home/ernie/analysis/cpptools -I/home/ernie/analysis/cpptools/CCfits -I/home/ernie/analysis/cpptools/cfitsio
INC = -I${APP3_CPP_INC} ${CPP_TOOLS_INC} -I${GSL_INC}

GSL_LIB = -L/home/ernie/analysis/cpptools/gsl-1.15/.libs/ -L/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs/
CCFITS_LIB = /home/ernie/analysis/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /home/ernie/analysis/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/home/ernie/analysis/cpptools/gsl-1.15/.libs/:/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs/:/home/ernie/analysis/cpptools/CCfits/.libs/:/home/ernie/analysis/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS = -Wall ${INC}
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS} -fopenmp

objs =  /home/ernie/analysis/apparatus3-analysis/qini/qini_utils.o /home/ernie/analysis/apparatus3-analysis/funcs/funcs.o /home/ernie/analysis/apparatus3-analysis/utils/utils.a 

testfits: testfits.o fits.a 
	g++  $^ ${objs} ${LFLAGS} -o testfits

fits.a: fits.h gaus1d.cpp gaus2d.cpp fermi2d.cpp fermiAzimuth.cpp
	indent fits.h
	indent gaus1d.cpp
	indent gaus2d.cpp
	indent fermi2d.cpp
	indent fermiAzimuth.cpp
	g++ ${CFLAGS} gaus1d.cpp gaus2d.cpp fermi2d.cpp fermiAzimuth.cpp -c 
	ar rcs fits.a gaus1d.o gaus2d.o fermi2d.o fermiAzimuth.o
	rm gaus1d.o gaus2d.o fermi2d.o fermiAzimuth.o

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




