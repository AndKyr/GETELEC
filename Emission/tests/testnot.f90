program testnot

use emission, only:  J_num_integ
use std_mat
implicit none


double precision, parameter:: F=1.2d0, W=4.5, R=100.d0, T=900.d0,gamma=10.d0
double precision::heat, Jcur
character:: regime


Jcur=J_num_integ(F,W,R,gamma,T,heat)

end program
