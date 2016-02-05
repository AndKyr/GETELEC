OBJECTS = emission.o std_mat.o
MODULES = emission.mod std_mat.mod
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none #-Wextra  -Wall -pedantic# -pedantic -O3
.PHONY: test jfplot

jfplot: test.exe pointsplot.py
	./test.exe
	python pointsplot.py

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
