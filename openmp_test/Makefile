all: openmp

INI_INC = /lab/software/apparatus3/cpp/qini
CCFITS_INC = /usr/local/include/CCfits
GSL_INC = /lab/software/apparatus3/cpp/gsl/include

INC = -I${INI_INC} -I${CCFITS_INC} -I${GSL_INC} 

GSL_LIB = /lab/software/apparatus3/cpp/gsl/lib

CFLAGS = -Wall ${INC} -fopenmp 
LFLAGS = -L${GSL_LIB} -lgsl -lgslcblas -lm -lCCfits -ltiff -fopenmp


objs =   


openmp: openmp.o 
	g++  $^ ${objs} ${LFLAGS} -o openmp

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




