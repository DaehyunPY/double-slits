# ver1 rk

TARGET = exe
OBJECTS = fftsg.o sde.o main.o 
COMMON_MOD = fftsg.f08 sde.f08 main.f08

F08 = /usr/local/bin/gfortran
FFLAGS = -fimplicit-none -fbounds-check

INCLUDE = -I/usr/local/include
LIB = -L/usr/local/lib
LAPACK = -llapack -lblas
FFW = -lfftw3 -lm





all: ${TARGET}

${TARGET}: ${OBJECTS}
	${F08} -o $@ ${OBJECTS}

${OBJECTS}: ${COMMON_MOD}
	${F08} -c ${COMMON_MOD} ${FFLAGS}

