program testnot

use emission, only:  J_num_integ
use std_mat
implicit none


double precision, parameter:: F=.d0, W=4.5d0, R=2.d0, T=600.d0,gamma=10.d0
double precision::heat, Jcur
character:: regime


Jcur=J_num_integ(F,W,R,gamma,T,heat)

end program