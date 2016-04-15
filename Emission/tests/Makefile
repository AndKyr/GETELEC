OBJECTS = std_mat.o pyplot_mod.o bspline.o
MODULES = std_mat.mod pyplot_mod.mod bspline.mod 
FFLAGS = -ffree-line-length-none #-fcheck=bounds#-Wextra -Wall -pedantic# -pedantic -O3

.PHONY: splinetest spline3d

spline3d: spline3d.exe
	./spline3d.exe

splinetest: splinetest.exe
	./splinetest.exe

%.exe: %.o $(OBJECTS)
	gfortran $<  $(FFLAGS) $(OBJECTS) -o $@

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90  $(MODULES) 
	gfortran $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.exe *.mod *.txt *.out *.pyf *~
