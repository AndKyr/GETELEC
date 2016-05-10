program surfacepoints

use interface_helmod, only: InterData, J_from_phi
use pyplot_mod , only: pyplot

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50

real(dp), allocatable   :: phi(:,:,:)
type(InterData)         :: this
integer                 :: icount, jcount, indsize(2),i

type(pyplot)            :: plt   !! pytplot handler


call read_phi(phi,this%grid_spacing)
this%grid_spacing = this%grid_spacing * 0.1d0


call J_from_phi(phi,this)

open(fidout,file='data/boundary_grid.xyz',action='write',status='replace')
write(fidout,*) size(this%Nstart,2)
write(fidout,*) 'eimaste treloi'

do i=1,size(this%Nstart,2)
    write(fidout,*) i, this%Nstart(:,i)*this%grid_spacing, this%Jem(i)
enddo 


close(fidout)

call plt%initialize(grid=.true.,xlabel='$1/F [nm/V]$',ylabel='$J (A/nm^2)$', &
            figsize=[20,15], font_size=font, title='FN-plot test', &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
call plt%add_plot(1.d0/this%F,log10(this%Jem),label='$current$', &
                    linestyle='b*',linewidth=2)
                    
call plt%add_plot(1.d0/this%F,log10(abs(this%heat)),label='$heat$', &
                    linestyle='r*',linewidth=2)



deallocate(this%Nstart, this%F, this%R,this%gamma, this%Jem, this%heat)

print * , 'Timing:'
print *, 'Cuurent:', this%TimeCur
print *, 'Set:', this%TimeInSet
print *, 'Fit', this%TimeFit
print *, 'Interpolate:', this%TimeInt 
    

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
    


