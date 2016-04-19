program surfacepoints

use emission, only: surf_points

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50

real(dp), allocatable   :: phi(:,:,:)
integer                 :: N, icount, jcount, indsize(2),i,j,k
integer, allocatable    :: inds(:,:)
real(dp), dimension(3)  :: grid_spacing,Ef

call read_phi(phi,grid_spacing)
grid_spacing=grid_spacing*0.1d0

inds=surf_points(phi)
indsize=shape(inds)
N=indsize(2)


open(fidout,file='data/boundary_grid.xyz',action='write',status='replace')
write(fidout,*) N
write(fidout,*) 'eimaste treloi'

do icount=1,N
    i=inds(1,icount)
    j=inds(2,icount)
    k=inds(3,icount)
    
    Ef(:)=([phi(i+1,j,k), phi(i,j+1,k), phi(i,j,k+1)] &
                    -[phi(i-1,j,k), phi(i,j-1,k), phi(i,j,k-1)])/grid_spacing
                
    if (norm2(Ef)>1.d0) then 
        jcount=jcount+1
    endif
    
    write(fidout,*) 0, i*grid_spacing(1), j*grid_spacing(2), &
             k*grid_spacing(3), norm2(Ef)
    
enddo

print *, 'jcount', jcount

close(fidout)

deallocate(phi,inds)
    

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
    


