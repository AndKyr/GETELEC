module std_mat
!this module contains some standard mathematical - matrix routines that
!are used quite frequently
implicit none

integer,parameter::dp=8

contains

subroutine csvprint(fileunit,dataarr)
	! write real numbers to a CSV file
	integer, intent(in):: fileunit
	double precision, intent(in) :: dataarr(:,:)
	integer            :: i, j 

	do i=1,size(dataarr,1)
		do j=1,size(dataarr,2)
			if (j==size(dataarr,2)) then
				write(fileunit,"(es22.15)",advance='no') dataarr(i,j)
			else
				write(fileunit,"(es22.15,',')",advance='no') dataarr(i,j)
			endif
		enddo
		write(fileunit,*) ''
	end do
end subroutine csvprint

subroutine csvread(fileunit,dataarr,rows,cols)
	! read real numbers from a CSV file
	integer, intent(in):: fileunit,rows,cols
	double precision, intent(out) :: dataarr(rows,cols)
	integer            :: i

	do i=1,rows
		read(fileunit,*) dataarr(i,:)
	end do
end subroutine csvread

pure function linspace(a,b,N) result(x)
	double precision,intent(in)::a,b
	integer,intent(in) ::N
	double precision :: dx,x(N)
	integer :: i
	
	dx=(b-a)/(N-1)
	do i=1,N
		x(i)=a+(i-1)*dx
	enddo
end function linspace

pure function logspace(a,b,N) result(x)
	double precision,intent(in)::a,b
	integer,intent(in) ::N
	double precision :: dlogx,logx(N),x(N)
	integer :: i
	
	dlogx=(log(b)-log(a))/(N-1)
	do i=1,N
		logx(i)=log(a)+(i-1)*dlogx
	enddo
	x=exp(logx)
end function logspace

pure function lininterp(yi,a,b,x) result(y)
!simple linear interpolation function
!appropriate for uniform linspace
	double precision, intent(in):: a,b,x, yi(:) ! yi interpolation array, a,b are x interval limits and x is the requested point
	integer :: Nnear,N
	double precision:: y,dx,dy,xnear
!	print *, yi
	
	if (x<a .or. x>b) then
		y=1.d200
		return
	endif
	N=size(yi)
	Nnear=nint((x-a)*N/(b-a))
	dx=(b-a)/dfloat(N-1)
	xnear=a+(Nnear-1)*dx
	if (x>xnear) then
		dy=yi(Nnear+1)-yi(Nnear)
	else
		dy=yi(Nnear)-yi(Nnear-1)
	endif
!	print*, 'Nnear=', Nnear, '|| N=', N, '|| dx=', dx, '|| dy=', dy, '|| xnear=', xnear, '|| x=', x
	y=yi(Nnear)+(dy/dx)*(x-xnear)
end function

pure function interp1(xi,yi,x) result(y)
!simple linear interpolation function (same as interp1 of matlab)
!appropriate for non-uniform linspace

!! be careful!! this  function is not thoroughly tested yet!!!!!
	double precision, intent(in) :: xi(:),yi(:),x
	double precision :: y
	integer :: i,j,Nlow
	
	do j=1,size(xi)
		if (xi(j)>x) then
			Nlow=j-1
			exit
		endif
	enddo
	if (Nlow<1 .or. j>size(xi)) then 
		y=1.d200
	else
		y=yi(Nlow)+((yi(Nlow+1)-yi(Nlow))/(xi(Nlow+1)-xi(Nlow)))*(x-xi(Nlow))
	endif
end function interp1

function diff2(f,x) result(y)
	double precision, intent(in)::x
	double precision, external::f
	double precision::y
	double precision,parameter::dx=1.d-2

	y=(f(x+dx)+f(x-dx)-2.d0*f(x))/(dx*dx)
end function diff2

  


   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

pure function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer, intent(in):: n
double precision, intent(in)::  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
end module std_mat

