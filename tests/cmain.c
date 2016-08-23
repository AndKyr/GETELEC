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
    data.mode = -20; data.full = 1;
    
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

int main(){
    
    int i, Nr = 32;
    double x[Nr], V[Nr], W = 4.5, Temp = 500.0, F = 5.0, R = 5.0, gamma = 10.0,
            Jem, heat;

    for(i=0; i<Nr; i++){
        x[i] = i*0.1+ 0.1;
        V[i] =  (F * R * x[i] * (gamma - 1.0) + 
                F * x[i] * x[i]) 
                / (gamma * x[i] + R * (gamma - 1.0));
    }
    
    libfun(W, Temp, Nr, x, V, &Jem, &heat);
    printf("\n I am C. \n Current: %4.7e, heat: %4.7e \n", Jem, heat);
    return 0;
}
