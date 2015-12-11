SHELL = /bin/sh

CC = gcc -std=c99 -fPIC
CFLAGS = -shared -dynamiclib -dy -fopenmp -Wall 
INCLUDES = -I. -I${SEISMIC_PATH}/src/API/c
LIBS =  -lm -lfftw3f

all: shotwem.dylib shotewem.dylib wavesep.dylib compute_angles.dylib

shotwem.dylib:shotwem.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} shotwem.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

shotewem.dylib:shotewem.c ewem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} shotewem.c ewem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

wavesep.dylib:wavesep.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} wavesep.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

compute_angles.dylib:compute_angles.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} compute_angles.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

clean:
	-rm -f *.o *.so *.dylib

