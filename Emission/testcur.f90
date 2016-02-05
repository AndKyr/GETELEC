program test


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: F=0.8d0, W=4.5d0, R=5.d0, T=1000.d0
double precision::J,G(4)
character :: regime





G=Gamow_general(F,W,R)
J=Cur_dens(F,W,R,T,regime)

print *, 'F=', F, '| Gam=', G(1), '|  dG/dE@Ef=' , G(2), '| dG/dE@Umax=', &
G(3), '|  Umax=', G(4) , '|  J= ', J , 'regime = ' , regime



end program
