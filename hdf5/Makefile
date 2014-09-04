MPICC=mpicc
MPIF90=mpif90
CC=gcc
FC=gfortran
CFLAGS=-g -I$(SCINET_HDF5_INC) -Wall

serial: 2darray 2darray-simple read2darray-simple f2darray-simple f2darray fread2darray-simple
parallel: parallel2darray fparallel2darray

HDF5CLIBS=-L$(SCINET_HDF5_LIB) -lhdf5 -L$(SCINET_SZIP_LIB) -lsz -lz
HDF5FLIBS=-L$(SCINET_HDF5_LIB) -lhdf5_fortran -lhdf5 -L$(SCINET_SZIP_LIB) -lsz -lz
HDF5INCS=-I$(SCINET_HDF5_INC) 

2darray: 2darray.o
	$(CC) -o $@ $< $(HDF5CLIBS) -lm

2darray-simple: 2darray-simple.o
	$(CC) -o $@ $< $(HDF5CLIBS) -lm

f2darray-simple: f2darray-simple.f90
	$(FC) -o $@ $< $(HDF5INCS) $(HDF5FLIBS) -lm

fread2darray-simple: fread2darray-simple.f90
	$(FC) -o $@ $< $(HDF5INCS) $(HDF5FLIBS) -lm

f2darray: f2darray.f90
	$(FC) -o $@ $< $(HDF5INCS) $(HDF5FLIBS) -lm

read2darray-simple: read2darray-simple.o
	$(CC) -o $@ $< $(HDF5CLIBS) -lm 

parallel2darray.o: parallel2darray.c
	$(MPICC) -g -c $< $(HDF5INCS)

parallel2darray: parallel2darray.o
	$(MPICC) -g -o $@ $< $(HDF5CLIBS) -lgpfs -lm

fparallel2darray.o: fparallel2darray.f90
	$(MPIF90) -g -c $< $(HDF5INCS)

fparallel2darray: fparallel2darray.o
	$(MPIF90) -g -o $@ $< $(HDF5FLIBS) -lgpfs -lm

clean: 
	rm -f *.o 2darray 2darray-simple read2darray-simple parallel2darray f2darray-simple f2darray fread2darray-simple fparallel2darray *.h5
