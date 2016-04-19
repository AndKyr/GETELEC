program splinetest

use bspline, only: db3ink,db3val
use pyplot_mod , only: pyplot
use std_mat, only: csvread, linspace

implicit none

integer, parameter      :: dp=8, fidin=8969, font=35

integer, parameter      :: kx=3, ky=3, kz=3, iknot=0, nline=64
real(dp), allocatable   :: fcn(:,:,:), tx(:), ty(:), tz(:)
integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
integer                 :: idx=0, idy=0, idz=0, iflag

integer                 :: nx, ny, nz, sz(3), i, j, k, istart, jstart, kstart
real(dp), allocatable   :: phi(:,:,:), x(:), y(:), z(:)
real(dp), dimension(3)  :: rstart, direc, Efstart, grid_spacing
real(dp), dimension(nline)  :: xline, yline, zline, rline, Vline

 character(len=2)       :: lnstyle(5)
 character(len=10)      :: labels


type(pyplot)            :: plt   !! pytplot handler
lnstyle=['r-', 'g-', 'b-', 'k-', 'm-']


call read_phi(phi,grid_spacing)


sz=shape(phi)
nx=sz(1)
ny=sz(2)
nz=sz(3)

call plt%initialize(grid=.false.,xlabel='r (nm)',ylabel='potential (V)', &
                figsize=[20,15], font_size=font, title='bspline test', &
                legend=.true.,axis_equal=.false., legend_fontsize=font, &
                xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)


allocate(x(nx),y(ny),z(nz),fcn(nx,ny,nz),tx(nx+kx),ty(ny+ky),tz(nz+kz))

x=[(grid_spacing(1)*(i-1), i=1,nx)]
y=[(grid_spacing(2)*(i-1), i=1,ny)]
z=[(grid_spacing(3)*(i-1), i=1,nz)]

rline=linspace(0.d0,30.d0,nline)
do j=-2,2
    istart=nx/2+j
    jstart=ny/2
    do kstart=1,nz
        if (phi(istart,jstart,kstart)>1.d-8) exit
    enddo
    
    print *, 'Nstart=',istart,jstart,kstart
    rstart=[x(istart),y(jstart),z(kstart-1)]
    Efstart=([phi(istart+1,jstart,kstart), &
             phi(istart,jstart+1,kstart), &
             phi(istart,jstart,kstart+1)]-phi(istart,jstart,kstart))/grid_spacing
    
    direc=Efstart/norm2(Efstart)

    print *, 'Direction=', direc
    !***************************************
    ! preparation of the data to interpolate
    call db3ink(x,nx,y,ny,z,nz,phi,kx,ky,kz,iknot,tx,ty,tz,fcn,iflag)
    if (iflag==0) print *, 'succesfully set parameters'
    !***************************************

    xline=rstart(1)+rline*direc(1)
    yline=rstart(2)+rline*direc(2)
    zline=rstart(3)+rline*direc(3)

    do i=1,nline
        call db3val(xline(i),yline(i),zline(i),idx,idy,idz,tx,ty,tz,nx,ny,nz, &
                    kx,ky,kz,fcn,Vline(i),iflag,inbvx,inbvy,inbvz,iloy,iloz)
    enddo


    call plt%add_plot(rline,Vline,label='$ Gamaw Sas $', &
                linestyle=lnstyle(j+3),linewidth=2)
        !*************************** present results
        if (iflag==0) print *, 'succesfully iterpolated'
enddo

deallocate(x,y,z,fcn,tx,ty,tz)
    
        
call plt%savefig('png/splinetest.png', pyfile='python/splinetest.py')


    !********************************

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

end program splinetest
    


