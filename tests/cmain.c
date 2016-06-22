#include <stdio.h>



extern int libfun(double W, double T, int Nr, double *xr, double *Vr, double *Jem,
            double *heat);

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
