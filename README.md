# GETELEC
General Tool for Electron Emission Calculations - A computational tool for calculating electron emission current and Nottingham effect heat for 

GETELEC is a scientific software for calculating electron emission currents and the Nottingham effect heat. For details see a 
forthcoming scientific publication.

Currently it is developed onlyfor Linux systems, but a version that can be compiled in other operating systems will be available in forthcoming
updates.

After downloading the zip file, the user needs to extract it and in the resulting folder execute ”make”
(GNU make is required) to compile it and build the GETELEC static and dynamic libraries. For successful building,
it is required that in the system there are installed gfortran (version 5 or later) and gcc (version 5 or later). By
executing ”make tests” some test programs are compiled and executed, outputting some of the plots included in the
present paper. For those tests to produce correct results, existing installation of python (version 2.7 or later) with
the accompanying libraries ”numpy” (version 0.13 or later), ”matplotlib” (version 1.3 or later) and ”scipy” (version
0.13 or later) is required.
