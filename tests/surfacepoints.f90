program surfacepoints

implicit none

integer, parameter      :: dp=8, fidout=8646, font=35, nline=50

real(dp), allocatable   :: phi(:,:,:), x(:), y(:), z(:), Ef(:,:)
integer                 :: sz(3), Nx, Ny, Nz, i, j, k, icount, N, jcount
integer, allocatable    :: ijk(:,:)
real(dp), dimension(3)  :: grid_spacing

real(dp)                ::F,R,gamma,W=4.5d0,T=1.d3,heat,Jcur
character               ::regime


call read_phi(phi,grid_spacing)
grid_spacing=grid_spacing*0.1d0

sz=shape(phi)
Nx=sz(1)
Ny=sz(2)
Nz=sz(3)
N=Nx*Ny*Nz
allocate(x(N), y(N), z(N), Ef(3,N), ijk(3,N))

icount=1
jcount=0
do k=2,Nz-1
    do j=2,Ny-1
        do i=2,Nx-1
            if (phi(i,j,k)<1.d-8 .and. ( &
                    phi(i+1,j,k)>1.d-8 .or. phi(i-1,j,k)>1.d-8 .or. &
                    phi(i,j+1,k)>1.d-8 .or. phi(i,j-1,k)>1.d-8 .or. &
                    phi(i,j,k+1)>1.d-8 .or. phi(i,j,k-1)>1.d-8)) &
                    then
                x(icount)=i*grid_spacing(1)
                y(icount)=j*grid_spacing(2)
                z(icount)=k*grid_spacing(3)
                Ef(:,icount)=([phi(i+1,j,k), phi(i,j+1,k), phi(i,j,k+1)] &
                            -[phi(i-1,j,k), phi(i,j-1,k), phi(i,j,k-1)])/grid_spacing
                
                if (norm2(Ef(:,icount))>1.d0) then 
                    jcount=jcount+1
                endif
                
                icount=icount+1
            endif
        enddo
    enddo
enddo
print *, 'jcount', jcount

open(fidout,file='data/boundary_grid.xyz',action='write',status='replace')
    
write(fidout,*) icount-1
write(fidout,*) 'eimaste treloi'
    
do i=1,icount-1
    write(fidout,*) 0,' ',x(i),' ',y(i),' ',z(i),' ',norm2(Ef(:,i))
enddo
    
close(fidout)

deallocate(phi,x,y,z,Ef,ijk)
    

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
    


