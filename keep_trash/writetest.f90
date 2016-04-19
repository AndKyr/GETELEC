program writetest

implicit none

integer :: i,fid=189
double precision :: x(100)

x=[(0.1d0*i, i=1,100)]

open(fid,file='skata.csv',action='read')!,status='replace')

read(fid,*) x

print *, x
close(fid)
end program writetest

