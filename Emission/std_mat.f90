module std_mat
!this module contains some standard mathematical - matrix routines that
!are used quite frequently
implicit none

integer,parameter::dp=8

contains

subroutine plot1(xdata,ydata)
	real(dp), intent(in):: xdata(:),ydata(:)

	integer, parameter:: fidx=2365689,fidy=24685945

	open(fidx,file='xdata.csv')
	open(fidy,file='ydata.csv')
	write(fidx,"(es22.15)") xdata
	write(fidy,"(es22.15)") ydata
	call system('python plot.py')
	close(fidx)
	close(fidy)
end subroutine

subroutine plot(xdata,ydata)
	real(dp), intent(in):: xdata(:,:),ydata(:,:)

	integer, parameter:: fidx=2365689,fidy=24685945

	open(fidx,file='xdata.csv')
	open(fidy,file='ydata.csv')
	call csvprint(fidx,xdata)
	call csvprint(fidy,ydata)
	call system('python ../plot.py')
	close(fidx)
	close(fidy)
end subroutine

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
	Nlow=0
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

function diff2(f,x) result(y)!second derivative of a function f at point x
	double precision, intent(in)::x
	double precision, external::f
	double precision::y
	double precision,parameter::dx=1.d-2

	y=(f(x+dx)+f(x-dx)-2.d0*f(x))/(dx*dx)
end function diff2

function local_min ( a, b, eps, t, f, x )
!*****************************************************************************80
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.  
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
!    perhaps non-global, minimum to the same accuracy.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the 
!    golden section step, 01 July 2013.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    01 July 2013
!
!  Author:
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real ( kind = 8 ) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real ( kind = 8 ) LOCAL_MIN, the value F(X).
!*****************************************************************************
  implicit none

  real(kind=8),intent(in):: a,b,eps,t
  real(kind=8),external:: f
  real(kind=8),intent(out)::x
  real(kind=8):: c,d,e,fu,fv,fw,fx,local_min,m,p,q,r,sa,sb,t2,tol,u,v,w
!
!  C is the square of the inverse of the golden ratio.
!
  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

  sa = a
  sb = b
  x = sa + c * ( b - a )
  w = x
  v = w
  e = 0.0D+00
  fx = f ( x )
  fw = fx
  fv = fw

  do
    m = 0.5D+00 * ( sa + sb ) 
    tol = eps * abs ( x ) + t
    t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
    if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
      exit
    end if
!
!  Fit a parabola.
!
    r = 0.0D+00
    q = r
    p = q

    if ( tol < abs ( e ) ) then

      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )

      if ( 0.0D+00 < q ) then
        p = - p
      end if

      q = abs ( q )

      r = e
      e = d

    end if

    if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
         q * ( sa - x ) < p .and. &
         p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
      d = p / q
      u = x + d
!
!  F must not be evaluated too close to A or B.
!
      if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

        if ( x < m ) then
          d = tol
        else
          d = - tol
        end if

      end if
!
!  A golden-section step.
!
    else

      if ( x < m ) then
        e = sb - x
      else
        e = sa - x
      end if

      d = c * e

    end if
!
!  F must not be evaluated too close to X.
!
    if ( tol <= abs ( d ) ) then
      u = x + d
    else if ( 0.0D+00 < d ) then
      u = x + tol
    else
      u = x - tol
    end if

    fu = f ( u )
!
!  Update A, B, V, W, and X.
!
    if ( fu <= fx ) then

      if ( u < x ) then
        sb = x
      else
        sa = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        sa = u
      else
        sb = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end do

  local_min = fx

  return
end

end module std_mat
