program surfacepoints

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50

real(dp), allocatable   :: phi(:,:,:), x(:), y(:), z(:)
integer                 :: sz(3), nx, ny, nz, i, j, k, icount
real(dp), dimension(3)  :: grid_spacing

real(dp)                ::F,R,gamma,W=4.5d0,T=1.d3,heat,Jcur
character               ::regime


call read_phi(phi,grid_spacing)
grid_spacing=grid_spacing*0.1d0

sz=shape(phi)
nx=sz(1)
ny=sz(2)
nz=sz(3)
allocate(x(nx),y(ny),z(nz))

icount=1
do i=1,nx
    do j=1,ny
        do k=1,nz
            if (phi(i,j,k)>1.d-8) exit
        enddo
        
        x(icount)=i*grid_spacing(1)
        y(icount)=j*grid_spacing(2)
        z(icount)=k*grid_spacing(3)
        icount=icount+1
    enddo
enddo

open(fidout,file='data/boundary_grid.xyz',action='write',status='replace')
    
write(fidout,*) size(x)
write(fidout,*) 'eimaste treloi'
    
do i=1,size(x)
    write(fidout,*) 0,' ',x(i),' ',y(i),' ',z(i), 1
enddo
    
close(fidout)

deallocate(phi,x,y,z)
    

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
    


