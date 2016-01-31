OBJECTS = emission.o
MODULES = emission.mod 
DEPS  = -llapack 
FFLAGS = -ffree-line-length-none #-Wextra   #-Wextra -pedantic -O3
.PHONY: test

test: test.exe
	./test.exe

%.exe: %.o $(OBJECTS)
	gfortran $< $(DEPS)  $(OBJECTS)  -o $@

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90  $(MODULES) 
	gfortran $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod *.csv *.txt
