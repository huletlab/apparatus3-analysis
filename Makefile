all: analyze  

GSL_INC = /home/ernie/analysis/cpptools/gsl-1.15/
APP3_CPP_INC = /home/ernie/analysis/apparatus3-analysis/
CPP_TOOLS_INC = -I/home/ernie/analysis/cpptools -I/home/ernie/analysis/cpptools/CCfits -I/home/ernie/analysis/cpptools/cfitsio
PYTHON_INC = /lab/software/epd-7.3-2-rh5-x86-64/include/python2.7/
#INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC} -I${PYTHON_INC}
INC = -I${APP3_CPP_INC} ${CPP_TOOLS_INC} -I${GSL_INC} 

GSL_LIB = -L/home/ernie/analysis/cpptools/gsl-1.15/.libs/ -L/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs/
CCFITS_LIB = /home/ernie/analysis/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /home/ernie/analysis/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/home/ernie/analysis/cpptools/gsl-1.15/.libs/:/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs/:/home/ernie/analysis/cpptools/CCfits/.libs/:/home/ernie/analysis/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS =  -Wall ${INC}
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS} -fopenmp

objs =   /home/ernie/analysis/apparatus3-analysis/funcs/funcs.o /home/ernie/analysis/apparatus3-analysis/qini/qini_utils.o /home/ernie/analysis/apparatus3-analysis/fits/fits.a /home/ernie/analysis/apparatus3-analysis/utils/utils.a
 
analyze: analyze.o ${objs}
	g++ analyze.o ${objs} ${LFLAGS} -o analyze
	chmod a+w analyze 
	cp -v analyze /home/ernie/analysis/apparatus3-analysis/bin/analyze
	rm analyze.o

analyze.o: Fermions.h
	indent Fermions.h
	indent analyze.cpp
	g++ ${CFLAGS} analyze.cpp -c
	
	

basler: basler.o ${objs} Fermions.h
	g++ basler.o ${objs} ${LFLAGS} -o basler
	chmod a+w basler
	cp -v basler /home/ernie/analysis/apparatus3-analysis/bin/basler

probe: probe.o ${objs} Fermions.h
	g++ probe.o ${objs} ${LFLAGS} -o probe
	chmod a+w probe
	cp -v probe /home/ernie/analysis/apparatus3-analysis/bin/probe

clean:
	rm -f *.o

.cpp.o: 
	indent $<
	g++ ${CFLAGS} $< -c




