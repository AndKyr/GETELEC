#include <stdio.h>

struct pass{
    int Nr;
    double *xr, *Vr;
};

extern void cur_dens_c(double *F, double *W, double *R, double *gamma, double *Temp,
                        int *mode, struct pass *, int *full,
                        double *Jem, double *heat, char *regime, int *sharp);
                        

int main(){
    
    double F, W, R, gamma, Temp, Jem, heat;
    int mode, full, sharp, i;
    char regime;
    struct pass arrays;
    
    F = 5.0; W = 4.5; R = 5.0; gamma = 10.0; Temp = 500.0; mode = -2; full = 1;
    
    arrays.Nr = 32;
    
    //arrays.xr = malloc(arrays.Nr * sizeof(double));
    //arrays.Vr = malloc(arrays.Nr * sizeof(double));
    
    double x[arrays.Nr], V[arrays.Nr]; 
    arrays.xr = x;
    arrays.Vr = V;  
    
    for(i=0; i<arrays.Nr; i++){
        arrays.xr[i] = i*0.1+ 0.1;
        arrays.Vr[i] =  (F * R * arrays.xr[i] * (gamma - 1.0) + 
                F * arrays.xr[i] * arrays.xr[i]) 
                / (gamma * arrays.xr[i] + R * (gamma - 1.0));
        //printf("%d\t%5.3e\t%5.3e\n", i, arrays.xr[i], arrays.Vr[i]);
    }
    
    cur_dens_c(&F, &W, &R, &gamma, &Temp, &mode, &arrays, &full, &Jem, &heat,
            &regime, &sharp);
            
    printf("\n I am C. \n Current: %4.7e, heat: %4.7e \n", Jem, heat);
    printf("Sharpness: %c, regime: %c\n", sharp, regime);
    //free(arrays.xr); free(arrays.Vr);
}
