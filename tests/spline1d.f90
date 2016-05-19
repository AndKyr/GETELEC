program splinetest

use bspline, only: db1ink,db1val
use pyplot_mod , only: pyplot
use std_mat, only: linspace

implicit none

integer, parameter 		:: dp=8
integer, parameter		:: nx=100,nxplot = 324, font = 32
real(dp), parameter		:: xmin=0.d0, xmax=6.d0, F=5.d0, R=5.d0, gamma = 10.d0
integer, parameter		:: kx=4, iknot=0

real(dp) 				:: tx(nx+kx), x(nx), y(nx), bcoef(nx), xplot(nxplot), yplot(nxplot)
integer 				:: inbvx, idx

integer					:: i,iflag


type(pyplot)            :: plt   !! pytplot handler

!**************initialize parameters not to be changed
inbvx = 	1
idx =   	0


!***************************************
! preparation of the data to interpolate
x=linspace(xmin,xmax,nx)
y = (F * R * x*(gamma - 1.d0) + F * x**2) &
                / (gamma * x + R * (gamma - 1.d0))

call db1ink(x,nx,y,kx,iknot,tx,bcoef,iflag)
if (iflag==0) print *, 'succesfully set parameters'
!***************************************

!***********************************
!create line of points to interpolate
xplot = linspace(xmin,xmax,nxplot)

call plt%initialize(grid=.true.,xlabel='r (nm)',ylabel='potential (V)', &
			figsize=[20,15], font_size=font, title='bspline test', &
			legend=.true.,axis_equal=.false., legend_fontsize=font, &
			xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)

do i=1,nxplot
	call db1val(xplot(i),idx,tx,nx,kx, bcoef,yplot(i),iflag,inbvx)
enddo
	!***********************************
call plt%add_plot(xplot,yplot,label='$interp$',linestyle='b-',linewidth=2)
call plt%add_plot(x,y,label='$data$',linestyle='r*',linewidth=2)

call plt%savefig('splinetest.png', pyfile='splinetest.py')


	!********************************

end program splinetest
	


