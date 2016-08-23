#include <stdio.h>
#include <stdlib.h>

extern int eval(const char *func, int nArgs, double **inReeal,
                 double **inImag, int blockSize, double *outReal, 
                double *outImag);

int main(){
    const int nArgs = 33, blockSize = 5;
    
    int i, j, iflag;
    double W[] = {3.5, 3.8, 4.2, 4.4, 4.5},
        Temp[] = {300., 400., 500., 600., 700.},
        F[]    = {3.5, 3.8, 4.2, 4.4, 4.5}, 
        R = 5.0, gamma = 10.0, x, V;
    
    double *inReal[nArgs], *inImag[nArgs];
    double *outReal, *outImag;
    
       
    for(i=0; i<nArgs; i++){
        inReal[i] = malloc(blockSize * sizeof(double));
        inImag[i] = malloc(blockSize * sizeof(double));
        for(j = 0; j<blockSize; j++){
            if (!i) {
                inReal[i][j] = Temp[i];
                inImag[i][j] = W[i];
            }
            else{
                x = i*0.1;
                inReal[i][j] = x;
                inImag[i][j] =  (F[j] * R * x * (gamma - 1.0) + F[j] * x * x) 
                        / (gamma * x + R * (gamma - 1.0));
            }
        }
    }
    
    outReal = malloc(blockSize * sizeof(double));
    outImag = malloc(blockSize  * sizeof(double));
        
    iflag = eval("getelec", nArgs, inReal, inImag, blockSize, outReal, outImag);
    for(i=0; i<blockSize; i++)
        printf("%d,\tCurrent: %4.7e, heat: %4.7e \n", i, outReal[i], outImag[i]);
    
    return 0;
}
