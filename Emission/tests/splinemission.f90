program splinetest

use pyplot_mod , only: pyplot
use std_mat, only: linspace
use emission, only: pot_interp, fitpot, Cur_dens

implicit none

integer, parameter      :: dp=8, fidin=8969, font=35, nline=50

real(dp), allocatable   :: phi(:,:,:), x(:), y(:), z(:)
integer                 :: Nstart(3), sz(3), nx,ny,nz,istart,jstart,kstart, j
real(dp), dimension(3)  :: grid_spacing
real(dp), dimension(nline)  :: rline, Vline, fitted

real(dp)                ::F,R,gamma,W=4.5d0,T=1.d3,heat,Jcur
character               ::regime

character(len=2)        :: lnstyle(5)
character(len=2)        :: mrkstyle(5)
character(len=10)       :: labels


type(pyplot)            :: plt   !! pytplot handler

mrkstyle=['r+', 'g+', 'b+', 'k+', 'm+']
lnstyle=['r-', 'g-', 'b-', 'k-', 'm-']


call read_phi(phi,grid_spacing)
grid_spacing=grid_spacing*0.1d0

call plt%initialize(grid=.false.,xlabel='r (nm)',ylabel='potential (V)', &
                figsize=[20,15], font_size=font, title='bspline test', &
                legend=.true.,axis_equal=.false., legend_fontsize=font, &
                xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)

sz=shape(phi)
nx=sz(1)
ny=sz(2)
nz=sz(3)

rline=linspace(0.d0,3.d0,nline)
do j=1,5
    istart=nx/2+j+2
    jstart=ny/2
    do kstart=1,nz
        if (phi(istart,jstart,kstart)>1.d-8) exit
    enddo
    
    Nstart=[istart,jstart,kstart]

    Vline=pot_interp(phi,grid_spacing,Nstart,rline)
    call fitpot(rline,Vline,F,R,gamma)
    
    fitted= (F*R*rline*(gamma-1.d0)+F*rline**2) / (gamma*rline+R*(gamma-1.d0))
    Jcur=Cur_dens(F,W,R,gamma,T,regime,heat)
    
    print *, 'F,R,gamma=', F,R,gamma
    print *, 'J,regime,heat=',Jcur,regime,heat

    call plt%add_plot(rline,Vline,label='$ Gamaw Sas $', &
                linestyle=mrkstyle(j),markersize=10)
                
    call plt%add_plot(rline,fitted,label='$ Gamaw Sas $', &
                linestyle=lnstyle(j),linewidth=1)
        !*************************** present results
enddo

deallocate(phi)
    
        
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
    


