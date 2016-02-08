
program test

use std_mat, only: lininterp,csvprint,linspace,csvread
implicit none


double precision, parameter:: xmaxVel=10.d0, F=5.d0,work=4.5d0, R=5.d0

integer, parameter :: NVel=200
double precision::Vel(NVel),x(Nvel), V(500),arrout(500,3),xnew(500)
integer:: i,fidout=1987

x=linspace(0.d0,xmaxVel,NVel)
Vel=F*R*x/(x+R)
xnew=linspace(0.d0,2.d0*R+5.d0,500)

do i=1,size(xnew)
	V(i)=ExtBar(xnew(i))
enddo


arrout(:,1)=xnew
arrout(:,2)=V
arrout(:,3)=0.d0

open(fidout,file="J-Fplot.csv",action="write",status="replace")
call csvprint(fidout,arrout)
close(fidout)


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










	
