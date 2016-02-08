
program test

use std_mat, only: lininterp,csvprint,linspace,csvread
implicit none


double precision, parameter:: xmaxVel=20.d0, F=5.d0,work=4.5d0, R=5.d0

integer, parameter :: NVel=201
double precision::Vel(NVel), Vappl,beta ,x(500),V(500),arrout(500,3)
integer:: i,fidout=1987,fidin=1953

open(fidin,file="potential.csv.in",action="read")
call csvread(fidin,Vel,201,1)
close(fidin)

beta=10.0*Vel(2)

Vappl=F/beta
Vel=Vappl*Vel
print *, lininterp(Vel,0.d0,xmaxVel,19.999999d0)

!x=linspace(0.d0,xmaxVel+5.d0,500)
!do i=1,size(x)
!	V(i)=ExtBar(x(i))
!enddo


!arrout(:,1)=x
!arrout(:,2)=V
!arrout(:,3)=0.d0

!open(fidout,file="J-Fplot.csv",action="write",status="replace")
!call csvprint(fidout,arrout)
!close(fidout)


contains

function ExtBar(x) result(V)!external barrier model
		double precision, intent(in) :: x
		double precision::V,dx,dVel,Velectr
		integer:: NVel
		if (x<xmaxVel) then
			Velectr=lininterp(Vel,0.d0,xmaxVel,x)
		else
			NVel=size(Vel)
			dx=xmaxVel/(NVel-1)
			dVel=Vel(NVel)-Vel(NVel-1)
			Velectr=Vel(NVel)+(dVel/dx)*(x-xmaxVel)
		endif
		V= work -Velectr -0.36d0/(x+(0.5d0*(x**2))/R)
end function ExtBar

end program










	
