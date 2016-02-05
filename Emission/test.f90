program test


use emission, only: Cur_dens,Gamow_general
use std_mat
implicit none


double precision, parameter:: Fmin=.5d0, Fmax=10.d0, W=4.5d0, &
R=10.d0, T=1000.d0
integer, parameter :: Nvals=100
double precision::J(Nvals),G(4),F(Nvals),dF,t1,t2,arrout(Nvals,3),regnum(Nvals)
character :: regime
integer:: i,fid=1987



F=1.d0/linspace(1.d0/Fmax,1.d0/Fmin, Nvals)
call cpu_time(t1)
do i=1,Nvals
	!G=Gamow_general(F,W,R)
	J(i)=Cur_dens(F(i),W,R,T,regime)
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

open(fid,file="J-Fplot.csv",action="write",status="replace")
call csvprint(fid,arrout)
close(fid)

print *, 'elapsed CPU time:', t2-t1, 'sec' 
end program
