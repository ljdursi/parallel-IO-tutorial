CC=gcc
F90=gfortran
FC=gfortran
MPICC=mpicc
MPIF90=mpif90
CFLAGS=-g -Wall -I$(SCINET_NETCDF_INC)
FFLAGS=-g -Wall -I$(SCINET_NETCDF_INC)
NETCDFCLIBS=-L$(SCINET_NETCDF_LIB) -lnetcdf -L$(SCINET_HDF5_LIB) -lhdf5_hl -lhdf5 -lgpfs -lz
NETCDFFLIBS=-I$(SCINET_NETCDF_INC) -L$(SCINET_NETCDF_LIB) -lnetcdff -lnetcdf -lnetcdff -L$(SCINET_HDF5_LIB) -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lgpfs -lz

serial: 2darray-simple read2darray-simple 2darray f2darray-simple f2darray fread2darray-simple
parallel: parallel2darray fparallel2darray

2darray: 2darray.o
	$(CC) -o $@ $< $(NETCDFCLIBS) -lm

2darray-simple: 2darray-simple.o
	$(CC) -o $@ $<  $(NETCDFCLIBS) -lm

f2darray-simple: f2darray-simple.f90
	$(F90) -o $@ $< $(FFLAGS) $(NETCDFFLIBS) -lm

fread2darray-simple: fread2darray-simple.f90
	$(F90) -o $@ $< $(FFLAGS) $(NETCDFFLIBS) -lm

f2darray: f2darray.f90
	$(F90) -o $@ $< $(FFLAGS) $(NETCDFFLIBS) -lm

read2darray-simple: read2darray-simple.o
	$(CC) -o $@ $<  $(NETCDFCLIBS) -lm

parallel2darray.o: parallel2darray.c
	$(MPICC) -g -c -Wall $< -I$(SCINET_NETCDF_INC) 

parallel2darray: parallel2darray.o
	$(MPICC) -g -o $@ $<  $(NETCDFCLIBS) -lm

fparallel2darray.o: fparallel2darray.f90
	$(MPIF90) -Wall -g -c $< -I$(SCINET_NETCDF_INC) 

fparallel2darray: fparallel2darray.o
	$(MPIF90) -g -o $@ $<  $(NETCDFFLIBS) -lm

clean: 
	rm -f *.o 2darray 2darray-simple read2darray-simple parallel2darray fread2darray-simple f2darray-simple f2darray fparallel2darray *.nc
