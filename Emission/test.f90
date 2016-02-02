program test


use emission, only: Gamow_general
implicit none


double precision::G(2),F=0.5d0
integer:: i

do i=1,50
	G=Gamow_general(F,4.5d0,7.d0)
	F=F+0.3d0
enddo


end program
