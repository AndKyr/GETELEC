
FC = gfortran-5
MODOBJ = emission.o
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -O3 #-Wall -pedantic# -pedantic -O3

.PHONY: main spectroscopy spline3d splinemission surfacepoints
.SECONDARY: MODOBJ


surfacepoints: surfacepoints.exe
	surfacepoints.exe
# 	ovito data/boundary_grid.xyz

%.exe: %.o $(MODOBJ)
	$(FC) $(FFLAGS) $^ $(DEPS) -o $@

%.o : $(MODOBJ) %.f90
	$(FC) $(FFLAGS) -c $(lastword $^) -o $@

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@ 

clean:
	rm -rf bin/* obj/* modules/obj/*
