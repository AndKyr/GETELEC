#ifndef GETELEC_H_
#define GETELEC_H_

#ifdef __cplusplus 
extern "C" {
#endif

struct emission{//struct for interoperability with fortran module getelec
    double F, W, R, gamma, Temp;//input parameters
    double Jem, heat; //ouptut parameters
    double *xr, *Vr;// input vectors
    char regime, sharp; //ouput characters showing regimes
    int Nr, approx, mode, ierr; //Nr: length of vectors xr ,Vt
                             //full: logical for full calculation
                             //mode, ierr: ints for defining mode and outputing error 
};

//external function from getelec fortran module doing all the connection
extern void c_wrapper(struct emission * , int );

//the three basic call subroutines of getelec working on struct emission
int cur_dens_c(struct emission *);
int print_data_c(struct emission *, int);
int plot_data_c(struct emission *);

//functions needed for comsol          
static const char *error = NULL; 
int init(const char *);
const char * getLastError(void);
int eval(const char *, int , const double **, const double **,
        int , double *, double *);//main function that comsol calls
        
#ifdef __cplusplus
}
#endif
        
#endif /* GETELEC_H_ */
