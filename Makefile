VERSION          = 1.1.9f4
BOPT		 = O

EXE =  Mix_picar3d-${VERSION}

all: $(EXE)

CPP		 = /lib/cpp
FC		 = ifort
CC       = icc 
COMMFLAGS        = -O3
#COMMFLAGS        = -O3 -ipo 
#COMMFLAGS        = -O0 -traceback -check 
FFLAGS = ${COMMFLAGS}

#LFLAGS          = /STACK:21737418240 
LFLAGS          = 
CPPFLAGS         =

FOBJS=$(wildcard IBM/*.o)

CFLAGS =-Wall -g
SOURCES =$(wildcard FEMcodes/*.cpp)
CPPOBJS =  $(SOURCES:.cpp=.o)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

	
${EXE}: $(FOBJS) $(CPPOBJS)
	${FC} -o ${EXE} $(FOBJS) $(CPPOBJS) -L/usr/lib -lstdc++


