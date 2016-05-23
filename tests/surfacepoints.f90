program surfacepoints

use new_interface, only: InterData, J_from_phi , interp_set, surf_points
use pyplot_mod , only: pyplot

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50


type(InterData)         :: this
integer                 :: icount, jcount, indsize(2),i
integer, allocatable    :: inds(:,:)
real(dp), allocatable   :: Jem(:), heat(:), F(:)

type(pyplot)            :: plt   !! pytplot handler

this%W = 4.5d0
this%kT = 0.07d0
call read_phi(this%phi,this%grid_spacing)
this%grid_spacing = this%grid_spacing * 0.1d0
call interp_set(this)
inds = surf_points(this%phi)
allocate(Jem(size(inds,2)), heat(size(inds,2)), F(size(inds,2)))

do i = 1,size(inds,2)
    this%Nstart = inds(:,i)
    call J_from_phi(this)
    Jem(i) = this%Jem
    heat(i) = this%heat
    F(i) =  norm2(this%F)
enddo

open(fidout,file='data/boundary_grid.xyz',action='write',status='replace')
write(fidout,*) size(Jem)
write(fidout,*) 'eimaste treloi'

do i=1,size(inds,2)
    write(fidout,*) i, inds(:,i)*this%grid_spacing, log10(Jem(i))
enddo 


close(fidout)

call plt%initialize(grid=.true.,xlabel='$1/F [nm/V]$',ylabel='$J (A/nm^2)$', &
            figsize=[20,15], font_size=font, title='FN-plot test', &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
call plt%add_plot(1.d0/F,log10(Jem),label='$current$', &
                    linestyle='b.',linewidth=2)
                    
call plt%add_plot(1.d0/F,log10(abs(heat)),label='$heat$', &
                    linestyle='r.',linewidth=2)
                    
call plt%savefig('png/surfacepoints.png', pyfile='python/surfacepoints.py')

deallocate(inds, F, heat, Jem, this%phi, this%bcoef, this%tx, this%ty, this%tz)

print * , 'Timing:'
print *, 'Set:', this%TimeInSet
print *, 'Interpolate:', this%TimeInt 
print *, 'Current:', this%TimeCur

contains

subroutine read_phi(phi,grid_spacing)

    implicit none

    double precision, allocatable, intent(out)  :: phi(:,:,:)
    double precision, intent(out)               :: grid_spacing(3)
    integer::i,j,k,sizes(3),Nx,Ny,Nz,fid=154967

    Nx=size(phi,1)
    Ny=size(phi,2)
    Nz=size(phi,3)
    
    open(fid,file='data/phi_grid.dat',action='read')
    
    read(fid,*) sizes
    read(fid,*) grid_spacing
    Nx=sizes(1)
    Ny=sizes(2)
    Nz=sizes(3)

    allocate(phi(Nx,Ny,Nz)) 
        
    do i=1,Nx
        do j=1,Ny
            read(fid,*) phi(i,j,:)
        enddo
        read (fid,*)
    enddo
    
    close(fid)
end subroutine read_phi


end program surfacepoints
    


