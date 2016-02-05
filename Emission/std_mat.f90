module std_mat

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

function linspace(a,b,N) result(x)
	double precision,intent(in)::a,b
	integer,intent(in) ::N
	double precision :: dx,x(N)
	integer :: i
	
	dx=(b-a)/(N-1)
	do i=1,N
		x(i)=a+(i-1)*dx
	enddo
end function linspace

function logspace(a,b,N) result(x)
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
end module std_mat

