program test


use emission, only: Cur_dens,Gamow_general
implicit none


double precision::J,F=0.7d0,t1,t2,G(4),W=4.5d0,R=20.d0,T=900.d0
character :: regime
integer:: i

call cpu_time(t1)

do i=1,50
	!G=Gamow_general(F,W,R)
	J=Cur_dens(F,W,R,T,regime)
	print *, 'F=',F,'|  J=', J ,' | regime: ', regime
	!print *, F,G
	F=F+0.3d0
enddo

call cpu_time(t2)
print *, 'elapsed CPU time:', t2-t1, 'sec' 
end program
