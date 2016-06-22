CC = gcc
FC = gfortran

MODOBJ = modules/obj/std_mat.o modules/obj/bspline.o \
  modules/obj/levenberg_marquardt.o modules/obj/pyplot_mod.o modules/obj/getelec.o \
  modules/obj/heating.o
  
DEPS  = -lslatec
FFLAGS = -ffree-line-length-none -fbounds-check -Imod -O3 -fPIC 
CFLAGS = -O3 -fPIC

LIBNAME=libgetelec.a
TEMPLIB = libtemp.a
LIBS = /usr/lib/libslatec.a
AR=ar rcs
LINKLIBS = ar -rcT

all: $(MODOBJ) #make into library
	$(AR) $(TEMPLIB) -o $(MODOBJ)
	$(LINKLIBS) $(LIBNAME) $(TEMPLIB) $(LIBS)
	
	
.PHONY: ccall current main spectroscopy
.SECONDARY: $(MODOBJ)



current: bin/current.exe
	./bin/current.exe 5.0 4.5 800.0 5.0 0 T

main: bin/main.exe
	./bin/main.exe

spectroscopy: bin/spectroscopy.exe
	./bin/spectroscopy.exe

lib/libgetelec.so:	cobj/c_interface.o lib/libgetelec.a 
	$(CC) -fPIC -shared -o $@ $^ -lslatec   
	
lib/libemission.a: lib/libgetelec.a
	$(LINKLIBS) $@ $< $(LIBS) 
	
lib/libgetelec.a: $(MODOBJ)
	$(AR) $@ -o $(MODOBJ)
	
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
