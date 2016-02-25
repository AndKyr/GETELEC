program test


use emission, only: Jcur
use std_mat
use omp_lib
implicit none


double precision, parameter:: Fmin=.1d0, Fmax=5.d0, W=4.5d0, &
R=3.d0, T=550.d0,gamma=10.d0
integer, parameter :: Nvals=8000
double precision::J(Nvals),F(Nvals),t1,t2,arrout(Nvals,3),regnum(Nvals), &
ww(Nvals),RR(Nvals),TT(Nvals),ggam(Nvals)
character :: regime(Nvals)
integer:: nthreads=1,i,fidout=1987

!$ nthreads = omp_get_num_procs()
!$ print * , 'threads: ', nthreads
!$ call omp_set_dynamic(.false.)
!$ call omp_set_num_threads(nthreads)


F=1.d0/linspace(1.d0/Fmax,1.d0/Fmin, Nvals)
ww=W
RR=R
TT=T
ggam=gamma

!print *, 'ww=', ww, 'RR=', RR, 'TT=', TT, 'ggam=',ggam
call cpu_time(t1)
!do i=1,Nvals
!!	Vel=F(i)*R*x/(x+R)
!	J(i)=Cur_dens(F(i),W,R,gamma,T,regime)
!	if (regime=='f') then
!		regnum(i)=-1.d0
!	elseif (regime=='t') then
!		regnum(i)=1.d0
!	else
!		regnum(i)=0.d0
!	endif
!	if (isnan(J(i))) then
!		print *, 'Nan at i=', i, 'F=', F(i)
!	endif
!enddo
J=Jcur(F,W,R,gamma,T,regime)

call cpu_time(t2)
!print *, 'J=', J
where (regime=='f') regnum=-1.d0
where (regime=='t') regnum=1.d0
where (regime=='i') regnum=0.d0
arrout(:,1)=1.d0/F
arrout(:,2)=log10(J)
arrout(:,3)=regnum

open(fidout,file="J-Fplot.csv",action="write",status="replace")
call csvprint(fidout,arrout)
close(fidout)

print *, 'elapsed CPU time:', t2-t1, 'sec'


contains



end program
