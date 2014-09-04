CC=mpicc
FC=mpif90
CFLAGS=-Wall -std=c99

EXES=helloworldf helloworldc helloworld-noncontigc helloworld-noncontigf writeatf writeatc writeatallc writeatallf fileviewf fileviewc

.SUFFIXES : .f90

% : %.f90
	$(FC) -o $@ $^ $(FFLAGS)

all: $(EXES)

sine: sine.c

sinef: sinef.f90

sineplot: sine
	mpirun -np 4 ./sine
	./dosineplot 

sinefplot: sinef
	mpirun -np 4 ./sinef
	./dosineplot 

clean:
	rm -rf *.o *.~ *.txt $(EXES) *.dat sinef sine
