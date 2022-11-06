#! /usr/bin/bash

gcc -c -fPIC -O3 integrator.c -o integrator.o
gfortran -c -fPIC -O3 dilog.f -o dilog.o

gcc -fPIC -shared integrator.o dilog.o -o libintegrator.so

rm *.o