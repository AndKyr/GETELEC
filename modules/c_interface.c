#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "getelec.h"

//external function from getelec fortran module doing all the connection
extern void c_wrapper(struct emission * data, int ifun);

//the three basic call subroutines of getelec working on struct emission
int cur_dens_c(struct emission *data){c_wrapper(data,0); return 0;}
int print_data_c(struct emission *data, int full){
    if (full) c_wrapper(data,2);//print in full mode
    else c_wrapper(data,1); //print only scalar data
    return 0;
}    
int plot_data_c(struct emission *data){c_wrapper(data,3); return 0;}


int print_C_data(struct emission *this){
    printf("F = %f\n", this->F);
    printf("W = %f\n", this->W);
    printf("R = %f\n", this->R);
    printf("gamma = %f\n", this->gamma);
    printf("Temp = %f\n", this->Temp);
    printf("Jem = %g\n", this->Jem);
    printf("heat = %g\n", this->heat);
    printf("regime = %d\n", this->regime);
    printf("sharp = %d\n", this->sharp);
    printf("Nr = %d\n", this->Nr);
    printf("approx = %d\n", this->approx);
    printf("mode = %d\n", this->mode);
    printf("ierr = %d\n", this->ierr);
    
    printf("   i\t   xr\t  Vr\n");
    
    for (int i = 0; i < this->Nr; ++i)
        printf("%d,  %g,  %g\n", i, this->xr[i], this->Vr[i]);

    return 0;
}

/***********************************************************************************
 * Comsol interface functions.
 * Interface for interoperability of GeTElEC with COMSOL. The comsol name of the
 * external function should be "getelec". 
 * nArgs-1 is the number of potential points given for each emission point. 
 * Each inReal array column contains the x points on the line outside each 
 * emission point (each emission point corresponds to a different column).
 *  
 * The first element of each column contains the temperature.
 * The same applies to inImag for the potential values. 
 * First element contains the 
 * work function of each point. 
 * 
 * The number of points is equal to blockSize and is 
 * the number of columns of the arrays.
 * 
 * On output the blockSize length vectors outReal and outImag contain the resulting
 * current density and nottingham heat correspondingly 
 * *********************************************************************************/           
int init(const char *str) { return 1; }
const char * getLastError() { return error; }
 
int eval(const char *func, int nArgs, const double **inReal, const double **inImag,
            int blockSize, double *outReal, double *outImag){
//Main calling function for comsol
    int i, j, iout;
    struct emission pass;
    double *x, *V;
    
    if (strcmp("getelec", func) != 0) {
        error = "Unknown function";
        printf("%s", error);
        return 0;
    }
    
    x = malloc((nArgs - 1) * sizeof(double)) ;
    V = malloc((nArgs - 1) * sizeof(double));

    for (i = 0; i < blockSize; i++) { //loop over blocksize different emission points
        pass.W = inImag[0][i]; //first imaginary Arg is the work function 
        pass.Temp = inReal[0][i]; //first real Arg is the temperature
        for(j = 1; j < nArgs; j++){
            x[j-1] = inReal[j][i];
            V[j-1] = inImag[j][i];
        }
        
        pass.xr = x;
        pass.Vr = V; 
        pass.mode = -21;
        pass.Nr = nArgs - 1;
        pass.approx = 1;

        cur_dens_c(&pass);
        
        outReal[i] = pass.Jem;
        outImag[i] = pass.heat;
        
    }
    free(x);
    free(V);    
    return 1;
}




