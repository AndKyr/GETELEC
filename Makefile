MODOBJ = modules/obj/std_mat.o modules/obj/bspline.o \
  modules/obj/levenberg_marquardt.o modules/obj/emission.o \
  modules/obj/pyplot_mod.o
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -fbounds-check -Imod #-Wall -pedantic# -pedantic -O3

.PHONY: main spectroscopy spline3d splinemission

splinemission: bin/splinemission.exe
	./bin/splinemission.exe

spline3d: bin/spline3d.exe
	./bin/spline3d.exe

spectroscopy: bin/spectroscopy.exe
	./bin/spectroscopy.exe

main: bin/main.exe
	./bin/main.exe

bin/%.exe: obj/%.o $(MODOBJ)
	gfortran $(FFLAGS) $^ $(DEPS) -o $@

obj/%.o : tests/%.f90
	gfortran $(FFLAGS) -c  $< -o $@

modules/obj/%.o : modules/%.f90
	gfortran $(FFLAGS) -Jmod -c $< -o $@ 

clean:
	rm -rf *~ *.o *.exe *.mod *.csv *.txt *.out *.pyf
