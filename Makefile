CC = gcc
FC = gfortran
AR=ar rcs
LINKLIBS = ar -rcT

MODOBJ = modules/obj/std_mat.o modules/obj/ellfuns_tabulated.o \
  modules/obj/bspline.o modules/obj/pyplot_mod.o modules/obj/getelec.o \
  modules/cobj/c_interface.o
  
DEPS  = -lslatec
FFLAGS = -fcheck=all -Imod -O3 -fPIC -Llib
CFLAGS = -fbounds-check -O3 -fPIC -Imodules #-Wall -Wextra

# Get list of object files, with paths
SLATEC_SRC := $(shell find lib/slatec/ -name '*.f')
SLATEC_OBJ := $(SLATEC_SRC:%.f=%.o)

PWD = $(shell pwd)

LIBSTATIC=lib/libgetelec.a
LIBDEPS = lib/libslatec.a
LIBSHARED = lib/dynamic/libgetelec.so

CINTERFACE = modules/cobj/c_interface.o
DIRS = bin cobj mod obj modules/obj modules/cobj png
	
.PHONY: tests varyingTemp ctest KXerror
.SECONDARY: *.o #$(MODOBJ)

all: $(DIRS) $(LIBSFULL) $(LIBSHARED)

tests: varyingTemp ctest KXerror plots fitIV  

$(DIRS):
	mkdir -p $(DIRS)
	

highFtest: bin/highFtest.exe
	./bin/highFtest.exe

varyingTemp: bin/varyingTemp.exe
	./bin/varyingTemp.exe	

errortest: bin/errortest.exe
	./bin/errortest.exe

ctest: bin/ctest.out
	./bin/ctest.out
	
KXerror: bin/KXerror.exe
	./bin/KXerror.exe
	
plots: bin/plots.exe
	./bin/plots.exe

current: bin/current.exe
	./bin/current.exe 5. 4. 1000.
	
fitIV: tests/fitIV.py
	./tests/fitIV.py
	
thetaSC: bin/thetaSC.exe
	./bin/theta.exe 1.e-6, 10, 500
	

$(LIBSHARED): $(CINTERFACE) $(LIBSTATIC) $(LIBDEPS)
	mkdir -p lib/dynamic/; $(FC) -fPIC -shared -o $@ $^   
	
$(LIBSTATIC): $(MODOBJ)
	$(AR) $@ $(MODOBJ)
	
lib/libslatec.a: $(SLATEC_OBJ)
	$(AR) $@ $(SLATEC_OBJ)
	
lib/slatec/obj/*.o: lib/slatec/src/*.f
	$(FC) $(FFLAGS) -c $< -o $@ 


bin/%.out: cobj/%.o $(LIBSHARED)
	$(CC) -L./lib -o $@ $^ #-lgetelec
	
bin/%.exe: obj/%.o $(MODOBJ)
	$(FC) $(FFLAGS) -Llib $^ $(DEPS) -o $@
	
cobj/%.o : tests/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ 

obj/%.o : $(MODOBJ) tests/%.f90
	$(FC) $(FFLAGS) -c $(lastword $^) -o $@

modules/obj/%.o : modules/%.f90
	$(FC) $(FFLAGS) -Jmod -c $< -o $@ 
	
modules/cobj/%.o : modules/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ 

clean:
	rm -rf bin/* obj/* cobj/* lib/libgetelec.so lib/libgetelec.a \
		lib/libemission.a modules/obj/* modules/cobj/* cobj/*.o libslat* 
