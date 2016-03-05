program testnot

use emission, only: Cur_dens
use std_mat
implicit none


double precision, parameter:: F=2.179300291545189d0, W=4.5, R=100.d0, T=900.d0,gamma=10.d0
double precision::heat, Jcur
character:: regime


Jcur=Cur_dens(F,W,R,gamma,T,regime,heat)
print *, 'F=', F, 'heat=', heat , 'regime=', regime, '  || Jcur=', Jcur
Jcur=Cur_dens(F,W,R,gamma,T,regime,heat)
print *, 'F=', F, 'heat=', heat , 'regime=', regime, '  || Jcur=', Jcur

end program
