all: testfuncs

GSL_INC = /home/ernie/analysis/cpptools/gsl-1.15/
APP3_CPP_INC = /home/ernie/analysis/apparatus3-analysis/
CPP_TOOLS_INC = /home/ernie/analysis/cpptools/
INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC}

GSL_LIB = -L/home/ernie/analysis/cpptools/gsl-1.15/.libs/ -L/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs

RUN_TIME_PATHS = -R/home/ernie/analysis/cpptools/gsl-1.15/.libs/:/home/ernie/analysis/cpptools/gsl-1.15/cblas/.libs/

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




