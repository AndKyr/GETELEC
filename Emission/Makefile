OBJECTS = std_mat.o emission.o
MODULES = std_mat.mod emission.mod 
DEPS  = -lslatec
FFLAGS = -O3 -ffree-line-length-none #-Wextra  #-Wall -pedantic# -pedantic -O3
.PHONY: test jfplot testgam testcur testinterp

testinterp: testinterp.exe
	./testinterp.exe

testcur: testcur.exe
	./testcur.exe

testgam: testgam.exe
	./testgam.exe

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
	rm -f *.o *.exe *.mod *.csv *.txt *.out *.pyf
