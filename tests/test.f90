program test


use emission, only: Cur_dens
use std_mat
implicit none


double precision, parameter:: Fmin=0.05d0, Fmax=7.d0, W=4.5d0, &
R=10.d0, T=700.d0,gamma=10.d0
integer, parameter :: Nvals=256
double precision   :: J(Nvals),F(Nvals),t1,t2,regnum(Nvals), &
ww(Nvals),RR(Nvals),TT(Nvals),ggam(Nvals),heat(Nvals)
character :: regime(Nvals)
integer:: i


F=linspace(Fmin,Fmax, Nvals)
ww=W
RR=R
TT=T
ggam=gamma


call cpu_time(t1)

do i=1,Nvals
	J(i)=Cur_dens(F(i),W,R,gamma,T,regime(i),heat(i))
	print *, J(i), heat(i), regime(i)
enddo

call cpu_time(t2)

! where (regime=='f') regnum=-1.d0
! where (regime=='t') regnum=1.d0
! where (regime=='i') regnum=0.d0
! 
! arrout(:,1)=1.d0/F
! arrout(:,2)=log10(J)
! arrout(:,3)=regnum
! arrout(:,4)=heat
! 
! open(fidout,file="J-Fplot.csv",action="write",status="replace")
! call csvprint(fidout,arrout)
! close(fidout)

print *, 'elapsed CPU time:', t2-t1, 'sec'




end program
