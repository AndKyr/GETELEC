program test


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: Fmin=.5d0, Fmax=10d0, W=4.5d0, &
R=10.d0, T=1200.d0,xmaxVel=2.d0*R
integer, parameter :: Nvals=1,NVel=1000
double precision::J(Nvals),F(Nvals),t1,t2,arrout(Nvals,3),regnum(Nvals),Vel(NVel),x(NVel)
character :: regime
integer:: i,fidout=1987


!
!9.995995995995995E-01

x=linspace(0.d0,xmaxVel,NVel)

F=1.d0/linspace(1.d0/Fmax,1.d0/Fmin, Nvals)

F=8.555279609488738d-1
call cpu_time(t1)
do i=1,Nvals
	Vel=F(i)*R*x/(x+R)
	J(i)=Cur_dens(F(i),W,R,T,regime)!,Vel,xmaxVel)
	if (regime=='f') then
		regnum(i)=-1.d0
	elseif (regime=='t') then
		regnum(i)=1.d0
	else
		regnum(i)=0.d0
	endif
enddo
call cpu_time(t2)
arrout(:,1)=F
arrout(:,2)=log(J)
arrout(:,3)=regnum

open(fidout,file="J-Fplot.csv",action="write",status="replace")
call csvprint(fidout,arrout)
close(fidout)

print *, 'elapsed CPU time:', t2-t1, 'sec'


contains



end program
