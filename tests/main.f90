program main

use emission, only: Jcur
use std_mat, only: csvread,csvprint
use omp_lib
implicit none

integer,parameter::dp=8,Nvals=500
real(dp):: params(4),F(Nvals),W,R,T,gamma,J(Nvals),regnum(Nvals),arrout(Nvals,3),heat(Nvals)
character :: regime(Nvals)
integer:: i,fidout=1987,fidin=1821,nthreads

open(fidin,file="paramin.csv",action="read")
read(fidin,*) params(:)
read(fidin,*) F(:)
close(fidin)

W=params(1)
R=params(2)
T=params(3)
gamma=params(4)

J=Jcur(F,W,R,gamma,T,regime,heat)

where (regime=='f') regnum=-1.d0
where (regime=='t') regnum=1.d0
where (regime=='i') regnum=0.d0
arrout(:,1)=1.d0/F
arrout(:,2)=log(J)
arrout(:,3)=regnum

open(fidout,file="J-Fplot.csv",action="write",status="replace")
call csvprint(fidout,arrout)
close(fidout)

end program
