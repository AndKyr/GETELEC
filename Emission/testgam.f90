program testgam


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: Fmin=.1d0, Fmax=2.d0, W=4.5d0, &
R=5.d0, T=1000.d0
integer, parameter :: Nvals=10
double precision::J(Nvals),G(4),F(Nvals)
character :: regime
integer:: i



F=linspace(Fmin,Fmax, Nvals)
do i=1,Nvals
	G=Gamow_general(F(i),W,R)
	print *, F(i), G
enddo

end program
