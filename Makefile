CC = gcc
FC = gfortran-5
AR=ar rcs
LINKLIBS = ar -rcT

MODOBJ = modules/obj/std_mat.o modules/obj/bspline.o \
  modules/obj/pyplot_mod.o modules/obj/getelec.o
  
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -fbounds-check -Imod -O3 -fPIC -Llib
CFLAGS = -O3 -fPIC -Imodules#-Wall -Wextra

LIBSTATIC=lib/libgetelec.a
LIBDEPS = lib/libslatec.a
LIBSFULL = lib/libemission.a
LIBSHARED = lib/libgetelec.so

CINTERFACE = modules/cobj/c_interface.o
DIRS = bin cobj mod obj modules/obj modules/cobj
	
.PHONY: ccall current main spectroscopy
.SECONDARY: *.o #$(MODOBJ)

all: $(DIRS) $(LIBSFULL) $(LIBSHARED)

$(DIRS):
	mkdir -p $(DIRS)

$(LIBSFULL): $(LIBSTATIC)
	$(LINKLIBS) $@ $< $(LIBDEPS)
	
varyingTemp: bin/varyingTemp.exe
	./bin/varyingTemp.exe	

fitFN: bin/fitFNplot.exe
	./bin/fitFNplot.exe

ctest: bin/ctest.out
	./bin/ctest.out
	
cmain: bin/cmain.out
	./bin/cmain.out

current: bin/current.exe
	./bin/current.exe 5.0 4.5 800.0 5.0 -21 T

main: bin/main.exe
	./bin/main.exe
	
KXerror: bin/KXerror.exe
	./bin/KXerror.exe
	

$(LIBSHARED): $(CINTERFACE) $(LIBSTATIC) 
	$(CC) -fPIC -shared -o $@ $^ $(DEPS)    
	
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
	
modules/cobj/%.o : modules/%.c
	$(CC) $(CFLAGS) $^ -c -o $@ 

clean:
	rm -rf bin/* obj/* lib/*.so lib/libgetelec.a lib/libemission.a modules/obj/* \
		cobj/*.o 
