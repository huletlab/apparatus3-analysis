all: analyze  

GSL_INC = /lab/software/apparatus3/cpptools/gsl-1.15/
APP3_CPP_INC = /lab/software/apparatus3/cpp/
CPP_TOOLS_INC = -I/lab/software/apparatus3/cpptools -I/lab/software/apparatus3/cpptools/CCfits -I/lab/software/apparatus3/cpptools/cfitsio
PYTHON_INC = /lab/software/epd-7.3-2-rh5-x86-64/include/python2.7/
#INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC} -I${PYTHON_INC}
INC = -I${APP3_CPP_INC} ${CPP_TOOLS_INC} -I${GSL_INC} 

GSL_LIB = -L/lab/software/apparatus3/cpptools/gsl-1.15/.libs/ -L/lab/software/apparatus3/cpptools/gsl-1.15/cblas/.libs/
CCFITS_LIB = /lab/software/apparatus3/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /lab/software/apparatus3/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/lab/software/apparatus3/cpptools/gsl-1.15/.libs/:/lab/software/apparatus3/cpptools/gsl-1.15/cblas/.libs/:/lab/software/apparatus3/cpptools/CCfits/.libs/:/lab/software/apparatus3/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS =  -Wall ${INC}
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS} -lgomp

objs =   /lab/software/apparatus3/cpp/funcs/funcs.o /lab/software/apparatus3/cpp/qini/qini_utils.o /lab/software/apparatus3/cpp/fits/fits.a /lab/software/apparatus3/cpp/utils/utils.a
 
analyze: analyze.o ${objs}
	g++ analyze.o ${objs} ${LFLAGS} -o analyze
	chmod a+w analyze 
	#cp -v analyze /lab/software/apparatus3/cpp/bin/analyze
	cp -v analyze /lab/software/apparatus3/cpp/bin/analyze_TF_err
	rm analyze.o

analyze.o: Fermions.h
	indent Fermions.h
	indent analyze.cpp
	g++ ${CFLAGS} analyze.cpp -c
	
	

basler: basler.o ${objs} Fermions.h
	g++ basler.o ${objs} ${LFLAGS} -o basler
	chmod a+w basler
	cp -v basler /lab/software/apparatus3/cpp/bin/basler

probe: probe.o ${objs} Fermions.h
	g++ probe.o ${objs} ${LFLAGS} -o probe
	chmod a+w probe
	cp -v probe /lab/software/apparatus3/cpp/bin/probe

clean:
	rm -f *.o

.cpp.o: 
	indent $<
	g++ ${CFLAGS} $< -c




