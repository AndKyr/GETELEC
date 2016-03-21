program main

use emission, only: J_num_integ, spectroscopy
use std_mat, only: csvread,csvprint
implicit none

integer,parameter::dp=8,fidin=1987
real(dp):: params(5),F,W,R,T,gamma,J,heat

spectroscopy=.true.
open(fidin,file="paramin.csv",action="read")
read(fidin,*) params(:)
close(fidin)

W=params(1)
R=params(2)
T=params(3)
gamma=params(4)
F=params(5)


J=J_num_integ(F,W,R,gamma,T,heat)



end program
