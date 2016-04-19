program testgam


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: Fmin=.345046d0, Fmax=0.34507d0, W=1.1869316394798384d0, &
R=10.d0, T=1000.d0
integer, parameter :: Nvals=1
double precision::J(Nvals),G(4),F,xi(4),yi(4),Um=-1.d20,xm=-1.d20
character :: regime
integer:: i

!yi=[2.d0, 7.d0, 6.d0, 7.d0]
!print *, lininterp(yi,0.d0,10.d0,1.d0)


!F=linspace(Fmin,Fmax, Nvals)
F=0.92637240356083084d0
print *, 'F=', F
G=Gamow_general(F,W,R,Um,xm)
!do i=1,Nvals
!	G=Gamow_general(F(i),W,R,Um,xm)
!	print *, F(i), G
!enddo

end program
