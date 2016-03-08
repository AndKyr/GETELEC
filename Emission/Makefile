OBJECTS = std_mat.o emission.o
MODULES = std_mat.mod emission.mod 
DEPS  = -lslatec
FFLAGS = -O3 -ffree-line-length-none#-Wextra  #-Wall -pedantic# -pedantic -O3
.PHONY: main

main: main.exe
	./main.exe

%.exe: %.o $(MODULES) $(OBJECTS)
	gfortran $<  $(FFLAGS) $(OBJECTS) -o $@ #$(DEPS) 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

%.o : %.f90  $(MODULES) 
	gfortran $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod *.csv *.txt *.out *.pyf
