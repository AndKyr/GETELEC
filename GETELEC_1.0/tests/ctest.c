#include <stdio.h>
#include <stdlib.h>
#include "getelec.h"
#include <string.h>

//extern int eval(const char *func, int nArgs, double **inReeal,
                 //double **inImag, int blockSize, double *outReal, 
                //double *outImag);

int main(){
    // const int nArgs = 33, blockSize = 5;
    
    // int i, j, iflag;
    // double W[] = {3.5, 3.8, 4.2, 4.4, 4.5},
    //     Temp[] = {300., 400., 500., 600., 700.},
    //     F[]    = {3.5, 3.8, 4.2, 4.4, 4.5}, 
    //     R = 5.0, gamma = 10.0, x, V;
    
    // double *inReal[nArgs], *inImag[nArgs];
    // double *outReal, *outImag;
    
       
    // for(i=0; i<nArgs; i++){
    //     inReal[i] = malloc(blockSize * sizeof(double));
    //     inImag[i] = malloc(blockSize * sizeof(double));
    //     for(j = 0; j<blockSize; j++){
    //         if (!i) {
    //             inReal[i][j] = Temp[i];
    //             inImag[i][j] = W[i];
    //         }
    //         else{
    //             x = i*0.1;
    //             inReal[i][j] = x;
    //             inImag[i][j] =  (F[j] * R * x * (gamma - 1.0) + F[j] * x * x) 
    //                     / (gamma * x + R * (gamma - 1.0));
    //         }
    //     }
    // }
    
    // outReal = malloc(blockSize * sizeof(double));
    // outImag = malloc(blockSize  * sizeof(double));
        
    // iflag = eval("getelec", nArgs, inReal, inImag, blockSize, outReal, outImag);
    // for(i=0; i<blockSize; i++){
    //     printf("Called comsol wrapper from C successfully. Results:\n");
    //     printf("%d,\tCurrent: %4.7e, heat: %4.7e \n", i, outReal[i], outImag[i]);
    // }
        
    // for(i=0; i<nArgs; i++){
    //     free(inReal[i]);
    //     free(inImag[i]);
    // }
    // free(outReal);
    // free(outImag);
    
    struct emission pass;
    pass.F = 10.;
    pass.W = 4.5;
    pass.R = 500.;
    pass.Temp = 300.;
    pass.gamma = 1.1;
    pass.mode = 0;
    pass.approx = 0;
    pass.voltage = 500.;
    // strcpy(pass.pfilename, "in/GetelecPar.in");

    pass.pfilename = (char*) malloc(strlen("in/GetelecPar.in")*sizeof(char));
    pass.pfile_length = strlen("in/GetelecPar.in");
    strncpy(pass.pfilename, "in/GetelecPar.in", pass.pfile_length);

    cur_dens_SC(&pass);
    
    print_data_c(&pass, 1);

    
    printf("theta_SC = %e\n", theta_SC(1.e-5, 1.e3, 10.));

    double G[10], Wmin, Wmax;

    export_gamow(5.,10.,20.,10, &Wmin, &Wmax, G);

    for (int i = 0; i < 10; i++)
        printf("%g ", G[i]);

    printf("\n%g %g \n", Wmin, Wmax);
    
    
    
    return 0;
}
