program testnot

use emission, only: Notting_gen
use std_mat
implicit none


double precision, parameter:: F=1.8199d0, W=4.5, R=100.d0, T=550.d0,gamma=10.d0
double precision::heat
character:: regime


heat=Notting_gen(F,W,R,gamma,T,regime)
print *, 'F=', F, 'heat=', heat , 'regime=', regime

end program
