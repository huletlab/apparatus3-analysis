all: qini wini

APP3_CPP_INC = /lab/software/apparatus3/cpp/
CPP_TOOLS_INC = /lab/software/apparatus3/cpptools/

INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC}

CFLAGS = -Wall ${INC}
LFLAGS = -lm

objs =  qini_utils.o

qini: ${objs} qini.cpp qini.o
	g++  qini.o ${objs} ${LFLAGS} -o qini 
	cp -v qini ../bin/

wini: ${objs} wini.cpp wini.o
	g++  wini.o ${objs} ${LFLAGS} -o wini 
	cp -v wini ../bin/

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS}  $< -c




