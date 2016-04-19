program splinetest

use bspline, only: db2ink,db2val
use pyplot_mod , only: pyplot
use std_mat, only: csvread, linspace

implicit none

integer, parameter 		:: dp=8, fidin=8969, font=35
integer, parameter		:: nx=100, ny=100, nline=100
real(dp), parameter		:: xmin=0.d0, xmax=6.d0, ymin=-2.d0, ymax=4.d0
real(dp), parameter		:: rmin=0.d0, rmax=3.d0, Rcirc=2.d0
integer, parameter		:: kx=4, ky=4, iknot=0

real(dp) 				:: tx(nx+kx),ty(ny+ky),fcn(nx,ny), theta(3)=[0.3d0,0.5d0,0.7d0]
integer 				:: inbvx,inbvy,iloy, idx, idy

real(dp) 				:: Vmesh(nx,ny), x(nx), y(ny), rdir(2)
real(dp) 				:: Vline(nline),xstart,ystart,ri(nline)
integer					:: i,j,iflag

character(len=2)		:: lnstyle(5)
character(len=10)		:: labels


type(pyplot)            :: plt   !! pytplot handler

!**************initialize parameters not to be changed
inbvx = 	1
inbvy = 	1
iloy =  	1
idx =   	0
idy =  		0

open(fidin,file="Vmesh.csv", status="old")

!***************************************
! preparation of the data to interpolate
x=linspace(xmin,xmax,nx)
y=linspace(ymin,ymax,ny)
call csvread(fidin,Vmesh,100,100)
call db2ink(x,nx,y,ny,Vmesh,kx,ky,iknot,tx,ty,fcn,iflag)
if (iflag==0) print *, 'succesfully set parameters'
!***************************************

!***********************************
!create line of points to interpolate
ri=linspace(rmin,rmax,nline)

call plt%initialize(grid=.true.,xlabel='r (nm)',ylabel='potential (V)', &
			figsize=[20,15], font_size=font, title='bspline test', &
			legend=.true.,axis_equal=.false., legend_fontsize=font, &
			xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)

lnstyle=['r-', 'g-', 'b-', 'k-', 'm-']

do j=1,size(theta)
	
	
	xstart=Rcirc*sin(theta(j))
	ystart = Rcirc*cos(theta(j))+ymin
	rdir=[sin(theta(j)), cos(theta(j))]
	x=xstart+ri*rdir(1)
	y=ystart+ri*rdir(2)

	write(labels,'(F4.2)') theta(j)
	print *, labels
	
	
	!**********************************

	!********************************* interpolate
	do i=1,nline
		call db2val(x(i),y(i),idx,idy,tx,ty,nx,ny,kx,ky,fcn,Vline(i),iflag,inbvx,inbvy,iloy)
	enddo
	!***********************************
	call plt%add_plot(ri,Vline,label='$\theta = '//trim(labels)//' $',linestyle=lnstyle(j),linewidth=2)
	!*************************** present results
	if (iflag==0) print *, 'succesfully iterpolated'
enddo
	
		
call plt%savefig('splinetest.png', pyfile='splinetest.py')


	!********************************

end program splinetest
	


