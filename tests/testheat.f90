program testheat

use heating, only: Potential, interp_set, get_heat, HeatData, heateq

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50


type(Potential)         :: poten
type(HeatData)          :: heat
integer                 :: i



call read_phi(poten%phi,poten%grid_spacing)
poten%grid_spacing = poten%grid_spacing * 0.1d0

allocate(heat%tempinit(size(poten%phi,3)), heat%hpower(size(poten%phi,3)), &
heat%tempfinal(size(poten%phi,3)))
heat%tempinit = 700.d0
call interp_set(poten)
call get_heat(heat,poten)

heat%Tbound = 700.d0
heat%maxtime = 1.d6
heat%dt = 1.d-1



call heateq(heat)

print *, 'temperature after ', heat%tempfinal


deallocate (heat%tempinit, heat%hpower, heat%tempfinal)
deallocate(poten%phi,poten%bcoef, poten%tx, poten%ty, poten%tz)

contains

subroutine read_phi(phi,grid_spacing)

    implicit none

    real(dp), allocatable, intent(out)  :: phi(:,:,:)
    real(dp), intent(out)               :: grid_spacing(3)
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


end program
    
