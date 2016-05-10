
FC = gfortran-5
MODOBJ = modules/obj/std_mat.o modules/obj/bspline.o \
  modules/obj/levenberg_marquardt.o modules/obj/emission.o \
  modules/obj/interface_helmod.o modules/obj/pyplot_mod.o
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -fbounds-check -Imod #-O3 #-Wall -pedantic# -pedantic -O3

.PHONY: main spectroscopy spline3d splinemission surfacepoints
.SECONDARY: MODOBJ

test: bin/test.exe
	./bin/test.exe

main: bin/main.exe
	./bin/main.exe

surfacepoints: bin/surfacepoints.exe
	./bin/surfacepoints.exe
# 	ovito data/boundary_grid.xyz

splinemission: bin/splinemission.exe
	./bin/splinemission.exe

spline3d: bin/spline3d.exe
	./bin/spline3d.exe

spectroscopy: bin/spectroscopy.exe
	./bin/spectroscopy.exe

bin/%.exe: obj/%.o $(MODOBJ)
	$(FC) $(FFLAGS) $^ $(DEPS) -o $@

obj/%.o : $(MODOBJ) tests/%.f90
	$(FC) $(FFLAGS) -c $(lastword $^) -o $@

modules/obj/%.o : modules/%.f90
	$(FC) $(FFLAGS) -Jmod -c $< -o $@ 

clean:
	rm -rf bin/* obj/* modules/obj/*
