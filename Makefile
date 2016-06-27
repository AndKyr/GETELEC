CC = gcc
FC = gfortran
AR=ar rcs
LINKLIBS = ar -rcT

MODOBJ = modules/obj/std_mat.o modules/obj/bspline.o \
  modules/obj/levenberg_marquardt.o modules/obj/pyplot_mod.o modules/obj/getelec.o \
  modules/obj/heating.o
  
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -fbounds-check -Imod -O3 -fPIC 
CFLAGS = -O3 -fPIC

LIBSTATIC=lib/libgetelec.a
LIBDEPS = /usr/lib/libslatec.a
LIBSFULL = lib/libemission.a
LIBSHARED = lib/libgetelec.so

CINTERFACE = cobj/inter_comsol.o


	
.PHONY: ccall current main spectroscopy
.SECONDARY: *.o #$(MODOBJ)

ctest: bin/ctest.out
	./bin/ctest.out

current: bin/current.exe
	./bin/current.exe 5.0 4.5 800.0 5.0 -21 T

main: bin/main.exe
	./bin/main.exe
	


spectroscopy: bin/spectroscopy.exe
	./bin/spectroscopy.exe

$(LIBSHARED): $(CINTERFACE) $(LIBSTATIC) 
	$(CC) -fPIC -shared -o $@ $^ $(DEPS)   
	
$(LIBSFULL): $(LIBSTATIC)
	$(LINKLIBS) $@ $< $(LIBDEPS) 
	
$(LIBSTATIC): $(MODOBJ)
	$(AR) $@ -o $(MODOBJ)
	
bin/%.out: cobj/%.o $(LIBSHARED)
	$(CC) -L./lib -Wl,-rpath=./lib -o $@ $^ #-lgetelec
	
bin/%.exe: obj/%.o $(MODOBJ)
	$(FC) $(FFLAGS) $^ $(DEPS) -o $@
	
cobj/%.o : tests/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ 

obj/%.o : $(MODOBJ) tests/%.f90
	$(FC) $(FFLAGS) -c $(lastword $^) -o $@

modules/obj/%.o : modules/%.f90
	$(FC) $(FFLAGS) -Jmod -c $< -o $@ 

clean:
	rm -rf bin/* obj/* lib/* modules/obj/* *.a *.so *.o
