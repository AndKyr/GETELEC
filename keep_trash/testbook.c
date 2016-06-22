#include <stdio.h>
#include <stdlib.h>

struct pass{
    int lenc, lenf;
    double *c, *f;
};

extern void foo(struct pass *arrays);

int main(){
    struct pass arrays;
    double cc[2], ff[3] ;
    
    arrays.lenc = 2 ;
    arrays.lenf = 3 ;
    
    //arrays.c = (double *) malloc(arrays.lenc);
    //arrays.f = (double *) malloc(arrays.lenf);
    arrays.c = cc;
    arrays.f = ff;
    
    arrays.c[0] = 1.0; arrays.c[1] = 2.0;
    arrays.f[0] = 1.0; arrays.f[1] = 2.0;
    
    foo(&arrays);
}
