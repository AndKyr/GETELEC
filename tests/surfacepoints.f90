program surfacepoints

use heating, only: Potential, PointEmission, HeatData, get_heat, &
                emit_atpoint, poten_create
use pyplot_mod , only: pyplot

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50


type(Potential)         :: poten
type(PointEmission)     :: point
type(HeatData)          :: heat

type(pyplot)            :: plt   !! pytplot handler

real(dp), allocatable   :: phi_in(:,:,:), z(:)
real(dp)                :: gs(3)
integer                 :: i

call read_phi(phi_in, gs)
call poten_create(poten, phi_in, gs * 0.1d0)

allocate( heat%tempinit(size(poten%phi, 3)), heat%tempfinal(size(poten%phi, 3)), &
            heat%hpower(size(poten%phi, 3)), heat%J_avg(size(poten%phi, 3)) )

heat%tempinit = 300.d0
heat%Tbound = 300.d0

call get_heat(heat,poten)
allocate(z(heat%tipbounds(1):heat%tipbounds(2)))
z = [(i * poten%grid_spacing(3), i=heat%tipbounds(1),heat%tipbounds(2))]

call plt%initialize(grid=.true.,xlabel='$h [nm]$',ylabel='$P_h [W/nm^3]$', &
            figsize=[20,15], font_size=font, title='heating', &
            legend=.true., axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font, ytick_labelsize=font, axes_labelsize=font)
            
call plt%add_plot(z,log10(heat%hpower(heat%tipbounds(1):heat%tipbounds(2))), &
            label='$heat$', linestyle='b.',linewidth=2)
                
                    
call plt%savefig('png/surfacepoints.png', pyfile='python/surfacepoints.py')


!print * , 'Timing:'
!print *, 'Set 3D:', this%timings(1), 's,  called', this%counts(1), 'times'
!print *, 'Interpolate 3D:', this%timings(2), 's,  called', this%counts(2), 'times'
!print *, 'Fitting:', this%timings(3), 's,  called', this%counts(3), 'times'
!print *, 'Set 1D:', this%timings(4), 's,  called', this%counts(4), 'times'
!print *, 'Current full model:', this%timings(5), 's, called', this%counts(5), 'times'
!print *, 'Current full interp:', this%timings(6), 's, called', this%counts(6), 'times'
!print *, 'Current GTF:', this%timings(7), 's,  called', this%counts(7), 'times'
!print *, 'Total time:', sum(this%timings), 's'


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
    


