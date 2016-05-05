program surfacepoints

use emission, only: surf_points, J_from_phi

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50

real(dp), allocatable   :: phi(:,:,:), Jcur(:), heat(:)
integer                 :: icount, jcount, indsize(2),i
integer, allocatable    :: inds(:,:)
real(dp), dimension(3)  :: grid_spacing,Ef
real(dp)    t1,t2, times(4)

call read_phi(phi,grid_spacing)
grid_spacing=grid_spacing*0.1d0


call cpu_time(t1)
inds=surf_points(phi)
call cpu_time(t2)
print *, 'Surface points:' , t2-t1

allocate(Jcur(size(inds,2)),heat(size(inds,2)))


open(fidout,file='data/boundary_grid.xyz',action='write',status='replace')
write(fidout,*) size(inds,2)
write(fidout,*) 'eimaste treloi'


Jcur=J_from_phi(phi,grid_spacing,inds,4.5d0,100.d0,heat,times)

do i=1,size(inds,2)
    write(fidout,*) i, inds(:,i)*grid_spacing, Jcur(i)
enddo 


close(fidout)

deallocate(phi,inds,Jcur,heat)
    

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
    


