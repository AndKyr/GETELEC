OBJECTS = std_mat.o pyplot_mod.o bspline.o
MODULES = std_mat.mod pyplot_mod.mod bspline.mod 
FFLAGS = -ffree-line-length-none #-Wextra -Wall -pedantic# -pedantic -O3

.PHONY: splinetest

splinetest: splinetest.exe
	./splinetest.exe

%.exe: %.o $(MODULES) $(OBJECTS)
	gfortran $<  $(FFLAGS) $(OBJECTS) -o $@

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90  $(MODULES) 
	gfortran $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.exe *.mod *.txt *.out *.pyf *~
