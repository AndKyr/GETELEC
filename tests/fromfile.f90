program fromfile

use GeTElEC, only: EmissionData, cur_dens, print_data, plot_barrier, dp
use std_mat, only: linspace

type(EmissionData)      :: this
integer                 :: i, Nx, fid = 878564
real(dp), allocatable   :: x(:), V(:)
    

open(fid, file = 'data/barrierdata.dat', action = 'read')

read(fid, *), Nx
allocate(x(Nx), V(Nx))

do i = 1,Nx
    read(fid,*), x(i), V(i)
enddo

this%xr = x
this%Vr = V
this%W = 4.5d0
this%kT = 0.025
this%mode = -20

!call plot_barrier(this)

call cur_dens(this)

call print_data(this)
call plot_barrier(this)

end program fromfile
