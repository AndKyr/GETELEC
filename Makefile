CC = gcc
FC = gfortran
MODOBJ = modules/obj/std_mat.o modules/obj/bspline.o \
  modules/obj/levenberg_marquardt.o modules/obj/pyplot_mod.o modules/obj/getelec.o \
  modules/obj/new_interface.o modules/obj/heating.o
  
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -fbounds-check -Imod -O3 #-Wall -pedantic# -pedantic -O3
CFLAGS = -Imod -O3

LIBNAME=libemission.a
TEMPLIB = libtemp.a
LIBS = /usr/lib/libslatec.a
AR=ar rcs
LINKLIBS = ar -rcT

all: $(MODOBJ) #make into library
	ar rcs $(LIBNAME) -o $(MODOBJ)
	
	
.PHONY: main spectroscopy spline3d splinemission surfacepoints current
.SECONDARY: $(MODOBJ)

ccall: bin/ccall.exe
	./bin/ccall.exe

testheat: bin/testheat.exe
	./bin/testheat.exe

current: bin/current.exe
	./bin/current.exe 5.0 4.5 800.0 5.0 0 T

testinterp: bin/testinterp.exe
	./bin/testinterp.exe

main: bin/main.exe
	./bin/main.exe

surfacepoints: bin/surfacepoints.exe
	./bin/surfacepoints.exe
# 	ovito data/boundary_grid.xyz

splinemission: bin/splinemission.exe
	./bin/splinemission.exe

spline1d: bin/spline1d.exe
	./bin/spline1d.exe

spectroscopy: bin/spectroscopy.exe
	./bin/spectroscopy.exe
	
bin/%.exe: obj/%.o $(MODOBJ)
	$(FC) $(FFLAGS) $^ $(DEPS) -o $@
	
obj/ccall.o : tests/c_interface.c
	$(CC) $(CFLAGS) $^ -c -o $@ 

obj/%.o : $(MODOBJ) tests/%.f90
	$(FC) $(FFLAGS) -c $(lastword $^) -o $@

modules/obj/%.o : modules/%.f90
	$(FC) $(FFLAGS) -Jmod -c $< -o $@ 

clean:
	rm -rf bin/* obj/* modules/obj/* *.a
