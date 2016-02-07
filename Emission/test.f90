program test


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: Fmin=2.d0, Fmax=10.d0, W=4.5d0, &
R=3.d0, T=300.d0
integer, parameter :: Nvals=2
double precision::J(Nvals),G(4),F(Nvals),dF,t1,t2,arrout(Nvals,3), &
regnum(Nvals),v(201),beta,Vappl(Nvals)
character :: regime
integer:: i,fidout=1987,fidin=1953

open(fidin,file="potential.csv",action="read")
call csvread(fidin,v,201,1)
close(fidin)

beta=10.0*v(2)

print *, beta

F=1.d0/linspace(1.d0/Fmax,1.d0/Fmin, Nvals)
Vappl=F/beta

print *, Vappl(1)*v(201)

call cpu_time(t1)
do i=1,Nvals
	!G=Gamow_general(F,W,R)
	J(i)=Cur_dens(F(i),W,R,T,regime)!,Vappl(i)*v,20.d0)
	if (regime=='f') then
		regnum(i)=-1.d0
	elseif (regime=='t') then
		regnum(i)=1.d0
	else
		regnum(i)=0.d0
	endif
enddo
call cpu_time(t2)
arrout(:,1)=1.d0/F
arrout(:,2)=log(J)
arrout(:,3)=regnum

open(fidout,file="J-Fplot.csv",action="write",status="replace")
call csvprint(fidout,arrout)
close(fidout)

print *, 'elapsed CPU time:', t2-t1, 'sec' 
end program
