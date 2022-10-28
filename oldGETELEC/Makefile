CC = gcc
FC = gfortran
AR=ar rcs
LINKLIBS = ar -rcT

MODOBJ = modules/obj/std_mat.o modules/obj/ellfuns_tabulated.o \
  modules/obj/bspline.o modules/obj/pyplot_mod.o modules/cobj/c_interface.o


FFLAGS = -fcheck=all -Imod -O3 -fPIC -Llib
CFLAGS = -fbounds-check -O3 -fPIC -Imodules #-Wall -Wextra

# Get list of object files, with paths
SLATEC_SRC := $(shell find lib/slatec/ -name '*.f')
SLATEC_OBJ := $(SLATEC_SRC:lib/slatec/src/%.f=lib/slatec/obj/%.o)

PWD = $(shell pwd)

LIBSTATIC=lib/libgetelec.a
LIBDEPS=lib/libslatec.a
LIBSHARED = lib/libgetelec.so

CINTERFACE = modules/cobj/c_interface.o
DIRS = bin cobj mod obj modules/obj lib lib/slatec/obj modules/cobj png

# Create version control for compiler
FORTRANERR = 
GNUVERSION := $(shell $(FC) -dumpversion)

ifeq ($(shell test $(GNUVERSION) -gt 9; echo $$?),0)
	FORTRANERR = -fallow-argument-mismatch
endif
	
.PHONY: tests varyingTemp ctest KXerror
.SECONDARY: *.o #$(MODOBJ)

all: $(DIRS) $(LIBSFULL) $(LIBSHARED) lib/libslatec.so python/libintegrator.so

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

python/%.f.o: python/%.f
	$(FC) -c -O3 -fPIC $^ -o $@

python/%.c.o: python/%.c
	$(CC) -c -O3 -fPIC $^ -o $@

python/libintegrator.so: python/integrator.c.o python/dilog.f.o
	$(CC) -fPIC -shared -o $@ $^
	

$(LIBSHARED): $(CINTERFACE) $(LIBSTATIC) $(LIBDEPS)
	$(FC) -fPIC -shared -o $@ $^
	
$(LIBSTATIC): obj/getelec.o lib/libslatec.a $(MODOBJ)
	$(AR) $@ $< $(MODOBJ)
	
lib/libslatec.a: $(SLATEC_OBJ)
	$(AR) $@ lib/slatec/obj/*.o

lib/libslatec.so: $(SLATEC_OBJ)
	$(FC) -fPIC -shared lib/slatec/obj/*.o -o $@
	
lib/slatec/obj/%.o: lib/slatec/src/%.f
	$(FC) $(FORTRANERR) -fPIC -O3 -w -c $< -o $@ 

bin/%.out: cobj/%.o $(LIBSHARED)
	$(CC) -L./lib -o $@ $^
	
bin/%.exe: obj/%.o obj/getelec.o $(LIBSHARED) $(MODOBJ)
	$(FC) $(FFLAGS) $^ $(DEPS) -o $@
	
cobj/%.o : tests/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ 
	
obj/getelec.o: modules/getelec.f90 $(MODOBJ)
	$(FC) -Jmod $(FFLAGS) -c $< -o $@

obj/%.o : $(MODOBJ) tests/%.f90
	$(FC) $(FFLAGS) -c $(lastword $^) -o $@

modules/obj/%.o : modules/%.f90
	$(FC) $(FFLAGS) -Jmod -c $< -o $@ 
	
modules/cobj/%.o : modules/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ 

clean:
	rm -rf bin/* obj/* cobj/* lib/libgetelec.so lib/libgetelec.a \
		modules/obj/* modules/cobj/* cobj/*.o libslat*

clean-all:
	rm -rf bin/* obj/* cobj/* lib/libgetelec.so lib/libgetelec.a \
		modules/obj/* modules/cobj/* cobj/*.o libslat* lib/libslat* \
		lib/slatec/obj/*
	
