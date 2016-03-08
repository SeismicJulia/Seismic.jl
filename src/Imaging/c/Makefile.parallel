SHELL = /bin/sh

CC = gcc -std=c99 -fPIC
CFLAGS = -shared -fopenmp -Wall 
INCLUDES = -I. -I/home/stanton/soft/fftw-3.3.4/include -I${SEISMIC_PATH}/src/API/c
LIBS =  -lm -L/home/stanton/soft/fftw-3.3.4/lib -lfftw3f

all: shotwem.so shotewem.so wavesep.so compute_angles.so

shotwem.so:shotwem.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} shotwem.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

shotewem.so:shotewem.c ewem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} shotewem.c ewem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

wavesep.so:wavesep.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} wavesep.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

compute_angles.so:compute_angles.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c
	${CC} ${INCLUDES} -o $@ ${CFLAGS} compute_angles.c wem.c ${SEISMIC_PATH}/src/API/c/seismic.c ${LIBS} 

clean:
	-rm -f *.o *.so

