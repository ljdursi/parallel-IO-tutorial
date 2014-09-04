CC=gcc
FC=gfortarn
MPICC=mpicc
MPIFC=mpif90

CFLAGS=-Wall
FFLAGS=-Wall 

override ADIOS_DIR := ${SCINET_ADIOS_BASE}
override ADIOS_FINC := $(shell ${ADIOS_DIR}/bin/adios_config -c -f) 
override ADIOS_CINC := $(shell ${ADIOS_DIR}/bin/adios_config -c) 
override ADIOS_FLIB := $(shell ${ADIOS_DIR}/bin/adios_config -l -f)
override ADIOS_CLIB := $(shell ${ADIOS_DIR}/bin/adios_config -l)

all: parallel2darray fparallel2darray

gwrite_ArrayData.ch: adios_global.xml
	gpp.py $<

gwrite_ArrayData.fh: fadios_global.xml
	gpp.py $<

parallel2darray : parallel2darray.c gwrite_ArrayData.ch
	${MPICC} -g ${ADIOS_INC} $< ${ADIOS_CLIB} -o $@ ${CFLAGS}

fparallel2darray : fparallel2darray.f90 gwrite_ArrayData.fh
	${MPIFC} -g ${ADIOS_INC} $< ${ADIOS_FLIB} -o $@ ${FFLAGS}

clean:
	rm -f *.ch *.o parallel2darray *h5 *nc *bp *fh
