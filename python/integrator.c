#include <math.h>
//[Work] + [kT] + + [self.Wmin] + [self.Wmax] + list(self.Gpoly) + list(self.dG) )

#define Energ xx[0]
#define Work xx[1]
#define kT xx[2]
#define Wmin xx[3]
#define Wmax xx[4]
#define dGmin xx[5]
#define dGmax xx[6]
#define poly(i) (xx[n - i - 1])
#define Npoly n - 7


double lFD(double E, double kboltzT){
    if (E > 20. * kboltzT)
        return exp(-E / kboltzT);
    else if(E < -20. * kboltzT)
        return -E / kboltzT;
    else
        return log(1. + exp(-E / kboltzT));
}


double Gfun(int n, double *xx){
    double G = poly(0);
    double W = Work - Energ;
    double Wlim;


    if (W > Wmax)
        Wlim = Wmax;
    else if (W < Wmin)
        Wlim = Wmin;
    else
        Wlim = W;
    
    double pow = Wlim;

    for (int i = 1; i < Npoly; i++){
        G += poly(i) * pow;
        pow *= Wlim;
    }

    if (W > Wmax)
        G += dGmax * (W - Wmax);
    else if (W < Wmin)
        G += dGmin * (W - Wmin);

    return G;
}

double intfun(int n, double *xx){

    double G = Gfun(n, xx);
    double D;

    if (G > 40.)
        D = exp(-G);
    else
        D = 1. / (1. + exp(G));

    return lFD(Energ, kT) * D;
}

double intfun_dbg(int n , double *xx){
    return intfun(n, xx);
}

double intfun_Pn(int n, double *xx){

    double G = Gfun(n, xx);
    double D;

    if (G > 40.)
        D = exp(-G);
    else
        D = 1. / (1. + exp(G));

    double expo = -exp(-Energ / kT);
    return D * (lFD(Energ, kT) * Energ - kT * dilog_(&expo)) ;
}



