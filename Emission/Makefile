OBJECTS = emission.o
MODULES = emission.mod 
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -O3 #-fPIC -Wextra  -Wall -pedantic# -pedantic -O3
.PHONY: test

test: test.exe
	./test.exe

%.exe: %.o $(MODULES) $(OBJECTS)
	gfortran $<  $(OBJECTS) $(DEPS) -o $@

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90  $(MODULES) 
	gfortran $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod *.csv *.txt *.out
