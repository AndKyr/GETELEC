module std_mat
!this module contains some standard mathematical - matrix routines that
!are used quite frequently

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
	double precision, intent(in):: a,b,x, yi(:) ! yi interpolation array, a,b are x interval limits and x is the requested point
	integer :: Nlow,N
	double precision:: y,dx,dy,xlow
	
	if (x<a .or. x>b) then
		y=1.d200
		return
	endif
	N=size(yi)
	Nlow=ceiling((x-a)*N/(b-a))
	dx=(b-a)/(N-1)
	dy=yi(Nlow+1)-yi(Nlow)
	xlow=dfloat((Nlow-1))/dfloat((N-1))

	y=yi(Nlow)+(dy/dx)*(x-xlow)
end function

pure function interp1(xi,yi,x) result(y)
!simple linear interpolation function (same as interp1 of matlab)
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
end module std_mat

