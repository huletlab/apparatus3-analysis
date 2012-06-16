all: testfuncs

GSL_INC = /lab/software/apparatus3/cpptools/gsl-1.15/
APP3_CPP_INC = /lab/software/apparatus3/cpp/
CPP_TOOLS_INC = /lab/software/apparatus3/cpptools/
INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC}

GSL_LIB = -L/lab/software/apparatus3/cpptools/gsl-1.15/.libs/ -L/lab/software/apparatus3/cpptools/gsl-1.15/cblas/.libs

RUN_TIME_PATHS = -R/lab/software/apparatus3/cpptools/gsl-1.15/.libs/:/lab/software/apparatus3/cpptools/gsl-1.15/cblas/.libs/

CFLAGS = -Wall ${INC} 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm -Xlinker ${RUN_TIME_PATHS} 


objs =  

testfuncs: testfuncs.o funcs.o
	g++  $^ ${objs} ${LFLAGS} -o testfuncs

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




