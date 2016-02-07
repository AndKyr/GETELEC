program testgam


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: Fmin=.1d0, Fmax=2.d0, W=4.5d0, &
R=5.d0, T=1000.d0
integer, parameter :: Nvals=10
double precision::J(Nvals),G(4),F(Nvals),xi(4),yi(4)
character :: regime
integer:: i

yi=[2.d0, 7.d0, 6.d0, 7.d0]
print *, lininterp(yi,0.d0,10.d0,1.d0)

!F=linspace(Fmin,Fmax, Nvals)
!do i=1,Nvals
!	G=Gamow_general(F(i),W,R)
!	!print *, F(i), G
!enddo

end program
