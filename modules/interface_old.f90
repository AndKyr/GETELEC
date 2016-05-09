module interface_helmod

implicit none

 contains

function surf_points(phi) result(inds2)

real(dp), intent(in)        :: phi(:,:,:)
integer , allocatable       :: inds2(:,:), inds(:,:)

integer                     :: i,j,k,Nx,Ny,Nz,N,icount, sz(3)

sz=shape(phi)
Nx=sz(1)
Ny=sz(2)
Nz=sz(3)
N=Nx*Ny*Nz
allocate(inds(3,N))

icount=1
do k=2,Nz-1
    do j=2,Ny-1
        do i=2,Nx-1
            if (phi(i,j,k)<1.d-8 .and. ( &
                    phi(i+1,j,k)>1.d-8 .or. phi(i-1,j,k)>1.d-8 .or. &
                    phi(i,j+1,k)>1.d-8 .or. phi(i,j-1,k)>1.d-8 .or. &
                    phi(i,j,k+1)>1.d-8 .or. phi(i,j,k-1)>1.d-8)) &
                    then
                inds(:,icount)=[i,j,k]
                
                icount=icount+1
            endif
        enddo
    enddo
enddo

allocate(inds2(3,icount-1))
inds2=inds(:,1:icount-1)
deallocate(inds)

end function surf_points

function J_from_phi(phi,grid_spacing,Nstart,W,T,heat,times) result(Jcur)
    use bspline, only: db3ink,db3val
    
    real(dp),intent(in)     :: phi(:,:,:),grid_spacing(3),T,W
    real(dp), intent(out)   :: heat(:),times(4)
    integer, intent(in)     :: Nstart(:,:)
    
    integer, parameter      :: kx=2, ky=2, kz=2, iknot=0, nline=64
    
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    real(dp), allocatable   :: fcn(:,:,:), tx(:), ty(:), tz(:), x(:), y(:), z(:)
    !the above are spline-related parameters
    
    integer                 :: nx, ny, nz, sz(3), i,j, istart, jstart, kstart, Npoints
    real(dp), dimension(3)  :: direc, Efstart
    real(dp), dimension(nline):: xline, yline, zline, V
    
    real(dp)                :: Jcur(size(Nstart,2)),ti(5)
    real(dp)                :: F,R,gamma, Vline(nline), rline(nline)
    real(dp)                :: tinterp=0.d0, tfit=0.d0, tcur=0.d0, t2,t1
    
    character               :: regime
    
    call cpu_time(t1)
    rline=linspace(0.d0,3.d0,nline)
    
    sz=shape(phi)
    nx=sz(1)
    ny=sz(2)
    nz=sz(3)
    Npoints=size(Nstart,2)
    

    allocate(x(nx),y(ny),z(nz),fcn(nx,ny,nz),tx(nx+kx),ty(ny+ky),tz(nz+kz))
    x=[(grid_spacing(1)*(i-1), i=1,nx)]
    y=[(grid_spacing(2)*(i-1), i=1,ny)]
    z=[(grid_spacing(3)*(i-1), i=1,nz)]
    call db3ink(x,nx,y,ny,z,nz,phi,kx,ky,kz,iknot,tx,ty,tz,fcn,iflag)
    call cpu_time(t2)
    print *, 'Interpolation set:', t2-t1
    do j=1,Npoints
    
        istart=Nstart(1,j)
        jstart=Nstart(2,j)
        kstart=Nstart(3,j)
    
    !find direction of the line (same as Efield direction)
        Efstart=([phi(istart+1,jstart,kstart), &
                 phi(istart,jstart+1,kstart), &
                 phi(istart,jstart,kstart+1)]-[phi(istart-1,jstart,kstart), &
                 phi(istart,jstart-1,kstart), &
                 phi(istart,jstart,kstart-1)])/grid_spacing
                 
        if (norm2(Efstart)>1.d0) then
            direc=Efstart/norm2(Efstart)
            !set the line of interpolation and interpolate
            xline=x(istart)+rline*direc(1)
            yline=y(jstart)+rline*direc(2)
            zline=z(kstart)+rline*direc(3)
            
            call cpu_time(t1)
            do i=1,nline
                call db3val(xline(i),yline(i),zline(i),idx,idy,idz,tx,ty,tz,nx,ny,nz, &
                        kx,ky,kz,fcn,Vline(i),iflag,inbvx,inbvy,inbvz,iloy,iloz)
            enddo
            call cpu_time(t2)
            tinterp=tinterp+t2-t1
            
            call cpu_time(t1)
            call fitpot(rline,Vline,F,R,gamma)
            call cpu_time(t2)
            tfit=tfit+t2-t1
            
            call cpu_time(t1)
!             print *, 'F,R,gamma=', F,R,gamma
            if (F>0.d0) then
                Jcur(j)=Cur_dens(F,W,R,gamma,T,regime,heat(j)) !calculate current
            else
                Jcur(j)=1.d-200
            endif
            call cpu_time(t2)
            tcur=tcur+t2-t1
        else
            Jcur(j)=1.d-200
        endif
        
    enddo

    print *, 'Interpolation: ', tinterp, 'sec'
    print *, 'Fitting:', tfit, 'sec'
    print *, 'Current:', tcur,'sec'

end function J_from_phi

function pot_interp(phi,grid_spacing,Nstart,r,timings) result(V)

    use bspline, only: db3ink,db3val
    
    real(dp), intent(in)    :: phi(:,:,:), grid_spacing(3), r(:)
    integer, intent(in)     :: Nstart(:)
    real(dp), intent(out)   :: timings(2)
    !everything comes in in nm
    
    integer, parameter      :: kx=2, ky=2, kz=2, iknot=0
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    real(dp), allocatable   :: fcn(:,:,:), tx(:), ty(:), tz(:), x(:), y(:), z(:)
    !the above are spline-related parameters
    
    integer                 :: nx, ny, nz, sz(3), i, nline, istart, jstart, kstart
    real(dp), dimension(3)  :: direc, Efstart
    real(dp), dimension(size(r))  :: xline, yline, zline, V
    
    real(dp)                :: t1,t2,t3
    
    call cpu_time(t1)
    !setup parameters and allocate memory for interpolating
    sz=shape(phi)
    nx=sz(1)
    ny=sz(2)
    nz=sz(3)
    nline=size(r)
    allocate(x(nx),y(ny),z(nz),fcn(nx,ny,nz),tx(nx+kx),ty(ny+ky),tz(nz+kz))
    x=[(grid_spacing(1)*(i-1), i=1,nx)]
    y=[(grid_spacing(2)*(i-1), i=1,ny)]
    z=[(grid_spacing(3)*(i-1), i=1,nz)]
    call db3ink(x,nx,y,ny,z,nz,phi,kx,ky,kz,iknot,tx,ty,tz,fcn,iflag)
    
    call cpu_time(t2)
    istart=Nstart(1)
    jstart=Nstart(2)
    kstart=Nstart(3)
    
    !find direction of the line (same as Efield direction)
    Efstart=([phi(istart+1,jstart,kstart), &
             phi(istart,jstart+1,kstart), &
             phi(istart,jstart,kstart+1)]-[phi(istart-1,jstart,kstart), &
             phi(istart,jstart-1,kstart), &
             phi(istart,jstart,kstart-1)])/grid_spacing
    direc=Efstart/norm2(Efstart)
!     print *, 'direc=' , direc
    !set the line of interpolation and interpolate
    xline=x(istart)+r*direc(1)
    yline=y(jstart)+r*direc(2)
    zline=z(kstart)+r*direc(3)
    
    do i=1,nline
        call db3val(xline(i),yline(i),zline(i),idx,idy,idz,tx,ty,tz,nx,ny,nz, &
                    kx,ky,kz,fcn,V(i),iflag,inbvx,inbvy,inbvz,iloy,iloz)
    enddo
    
    call cpu_time(t3)
    
    timings=[t2-t1,t3-t2]
    
    deallocate(x,y,z,fcn,tx,ty,tz)
    
end function pot_interp

subroutine fitpot(x,V,F,R,gamma)
    !Fits V(x) data to F,R,gamma standard potential using L-V
    !minimization module
    use Levenberg_Marquardt, only: nlinfit
    
    real(dp), intent(in)    ::x(:),V(:)
    real(dp), intent(out)   ::F,R,gamma
    real(dp)                ::p(3),F2,Fend,var
    
    integer                 :: Nstart=4

    
    p(1) = (V(Nstart)-V(Nstart-1))/(x(Nstart)-x(Nstart-1))
    F2 = (V(Nstart+1)-V(Nstart))/(x(Nstart+1)-x(Nstart))
    Fend = (V(size(V))-V(size(V)-1))/(x(size(x))-x(size(x)-1))
    p(2) = abs(2.d0/((F2-p(1))/(x(Nstart)-x(Nstart-1))))
    p(3) = p(1)/Fend
!     print *, 'F,R,gamma first guess =', p
    var=nlinfit(fun,x,V,p)
    F=p(1)
    R=p(2)
    gamma=p(3)

    contains

    pure function fun(x,p) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: p(:)
        real(dp)             :: y
        y=(p(1)*p(2)*x*(p(3)-1.d0)+p(1)*x**2) / (p(3)*x+p(2)*(p(3)-1.d0))
    end function fun
end subroutine fitpot