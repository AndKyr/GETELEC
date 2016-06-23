/* A C interface for interoperability of GeTElEC with COMSOL. The comsol name of the
 * external function should be "getelec". 
 * 
 * nArgs-1 is the number of potential points
 * given for each emission point. Each inReal array column contains the x points
 * on the line outside each emission point (each emission point corresponds to 
 * a different column). The first element of each column contains the temperature.
 * The same applies to inImag for the potential values. First element contains the 
 * work function of each point. The number of points is equal to blockSize and is 
 * the number of columns of the arrays.
 * 
 * On output the blockSize length vectors outReal and outImag contain the resulting
 * current density and nottingham heat correspondingly*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


struct emission{
    double F, W, R, gamma, Temp;//input parameters
    double Jem, heat; //ouptut parameters
    double *xr, *Vr;// input vectors
    char regime, sharp; //ouput characters showing regimes
    int Nr, full, mode; //length of vectors xr ,Vt, and logical for full calculation 
};

extern void cur_dens_c(struct emission * data);
              
static const char *error = NULL;

int init(const char *str) { return 1; }

const char * getLastError() { return error; }

int eval(const char *func, int nArgs, const double **inReal, const double **inImag,
            int blockSize, double *outReal, double *outImag)
{
    int i, j;
    struct emission pass;
    double x[nArgs - 1], V[nArgs - 1];
    
    if (strcmp("getelec", func) != 0) {
        error = "Unknown function";
        printf("%s", error);
        return 0;
    }

    for (i = 0; i < blockSize; i++) { //loop over blocksize different emission points
        pass.W = inImag[0][i]; //first imaginary Arg is the work function 
        pass.Temp = inReal[0][i]; //first real Arg is the temperature
        for(j = 1; j < nArgs; j++){
            x[j-1] = inReal[j][i];
            V[j-1] = inImag[j][i];
        }
        pass.xr = x;
        pass.Vr = V; 
        pass.mode = -2;
        pass.Nr = nArgs - 1;
        pass.full = 1;
        
        cur_dens_c(&pass);
        
        outReal[i] = pass.Jem;
        outImag[i] = pass.heat;
    }
    
    return 1;
}
