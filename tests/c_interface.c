#include <stdio.h>

struct emission{
    double F, W, R, gamma, Temp;//input parameters
    double Jem, heat; //ouptut parameters
    double *xr, *Vr;// input vectors
    char regime, sharp; //ouput characters showing regimes
    int Nr, full, mode; //length of vectors xr ,Vt, and logical for full calculation 
};

extern void cur_dens_c(struct emission *data);

int libfun(double W, double T, int Nr, double *x, double *V,
            double *Jem, double *heat){
    
    struct emission data;
    int i;
    
    
    data.F = 5.0; data.R = 5.0; data.gamma = 10.0;
    data.W = W;   data.Temp = T;
    data.mode = -2; data.full = 1;
    
    data.Nr = Nr;
    
    //arrays.xr = malloc(arrays.Nr * sizeof(double));
    //arrays.Vr = malloc(arrays.Nr * sizeof(double));
    data.xr = x;
    data.Vr = V;
    

    cur_dens_c(&data);
    *Jem = data.Jem;
    *heat = data.heat;
    return 0;
}



