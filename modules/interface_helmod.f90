module interface_helmod

implicit none

integer, parameter  :: dp=8, Nr=32
real(dp), parameter :: Jlim=1.d-22, Fmin = 0.5d0

type, public       :: InterData

    real(dp), allocatable   :: F(:), R(:), gamma(:), Jem(:), heat(:)
                            !list of parameters defining barrier
    real(dp)                :: W=4.5d0, kT=0.05d0, grid_spacing(3) ! General external parameters
    integer, allocatable    :: Nstart(:,:) 
                            !indices of starting surface surf_points
    real(dp)                :: rline(Nr), grid(3) 
                            !length of V lines and grid spacing
    real(dp)                :: TimeCur=0.d0, TimeInSet=0.d0, TimeFit=0.d0, TimeInt=0.d0 
                            !timing variables
end type InterData
    

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

subroutine J_from_phi(phi,this)
    use bspline, only: db3ink,db3val
    use std_mat, only: linspace
    
    type(InterData), intent(inout)  :: this
    real(dp), intent(in)            :: phi(:,:,:)
    
    integer, parameter      :: kx=2, ky=2, kz=2, iknot=0
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    real(dp), allocatable   :: fcn(:,:,:), tx(:), ty(:), tz(:), x(:), y(:), z(:)
    !the above are spline-related parameters
    
    integer                 :: nx,ny,nz, i, j, istart,jstart,kstart, Npoints
    real(dp), dimension(3)  :: direc, Efstart
    real(dp), dimension(Nr) :: xline, yline, zline, Vline
    
    logical                 :: plot = .false.
    
    real(dp)                :: t1, t2

    
    call cpu_time(t1)
    this%rline = linspace(0.d0,2.d0,Nr)
    
    nx = size(phi,1)
    ny = size(phi,2)
    nz = size(phi,3)
    
    this%Nstart = surf_points(phi)
    Npoints = size(this%Nstart,2)
    allocate(this%F(Npoints),this%R(Npoints),this%gamma(Npoints), &
            this%Jem(Npoints), this%heat(Npoints))
    

    allocate(x(nx),y(ny),z(nz),fcn(nx,ny,nz),tx(nx+kx),ty(ny+ky),tz(nz+kz))
    x = [(this%grid_spacing(1)*(i-1), i=1,nx)]
    y = [(this%grid_spacing(2)*(i-1), i=1,ny)]
    z = [(this%grid_spacing(3)*(i-1), i=1,nz)]
    call db3ink(x,nx,y,ny,z,nz,phi,kx,ky,kz,iknot,tx,ty,tz,fcn,iflag)
    call cpu_time(t2)
    this%TimeInSet = t2 - t1
    
    do j=1,Npoints
        istart = this%Nstart(1,j)
        jstart = this%Nstart(2,j)
        kstart = this%Nstart(3,j)
    
        !find direction of the line (same as Efield direction)
        ! always 
        Efstart=([phi(istart+1,jstart,kstart), &
                 phi(istart,jstart+1,kstart), &
                 phi(istart,jstart,kstart+1)]-[phi(istart-1,jstart,kstart), &
                 phi(istart,jstart-1,kstart), &
                 phi(istart,jstart,kstart-1)])/this%grid_spacing
                 
        if (norm2(Efstart)> Fmin) then
                 
            direc=Efstart/norm2(Efstart)
                !set the line of interpolation and interpolate
            xline=x(istart)+this%rline*direc(1)
            yline=y(jstart)+this%rline*direc(2)
            zline=z(kstart)+this%rline*direc(3)
                
            call cpu_time(t1)
            do i=1,Nr
                call db3val(xline(i),yline(i),zline(i),idx,idy,idz,tx,ty,tz,nx,ny,nz, &
                        kx,ky,kz,fcn,Vline(i),iflag,inbvx,inbvy,inbvz,iloy,iloz)
            enddo
            call cpu_time(t2)
            this%TimeInt = this%TimeInt + t2 - t1
            
!            if (j==4196) then
!                plot = .true.
!            else
!                plot = .false.
!            endif
                
!            call cpu_time(t1)
!            call fitpot(this%rline,Vline,this%F(j),this%R(j),this%gamma(j),plot)
!            call cpu_time(t2)
!            this%TimeFit = this%TimeFit + t2 - t1
        else
            this%F(j) = norm2(Efstart)
            this%R(j) = 1.d4
            this%gamma(j) = 1.d0
        endif
            
        call cpu_time(t1)
        call current(this,j)
        call cpu_time(t2)
        this%TimeCur = this%TimeCur + t2 - t1
        
    enddo
end subroutine J_from_phi

subroutine current(this,i)
!calls emission module and calculates emission for the ith 
    use emission, only: EmissionData, cur_dens, print_data
    
    type(InterData), intent(inout)  :: this
    integer, intent(in)             :: i
    
    type(EmissionData)               :: emit

    emit%F = this%F(i)
    emit%W = this%W
    emit%R = abs(this%R(i))
    emit%gamma = this%gamma(i)
    emit%kT = this%kT
    emit%full = .false.
    
    call cur_dens(emit)
    
    if (emit%Jem > Jlim) then
        emit%full = .true.
        call cur_dens(emit)
    endif
    
    if (isnan(emit%Jem) .or. abs(emit%Jem) > 1.d201) then
        call print_data(emit)
    endif
    
    this%Jem(i) = emit%Jem
    this%heat(i) = emit%heat
    
end subroutine current

!subroutine fitpot(x,V,F,R,gamma,plot) 
!!Fits V(x) data to F,R,gamma standard potential using nr fitting module
!    use pyplot_mod, only: pyplot
!    use fitting, only: mrqmin

!    integer, parameter      :: Nstart = 3, Np = 3
!    real(dp), intent(in)    :: x(:), V(:)
!    logical, intent(in), optional   :: plot
!    real(dp), intent(out)   :: F,R,gamma
    
    
!    real(dp)                :: sig(size(x)), Vfitted(size(x)), dfit(size(x),3)
!    real(dp)                :: a(Np), F2, Fend, chisq, alamda = -0.1d0
!    real(dp)                :: covar(Np,Np), alpha(Np,Np)
!    logical                 :: mask(Np) = [.true., .true., .true.]
    
!    type(pyplot)            :: plt
  
!    sig(1:) = 1.d-1
    
    
!    a(1) = (V(Nstart)-V(Nstart-1))/(x(Nstart)-x(Nstart-1))
!    F2 = (V(Nstart+1)-V(Nstart))/(x(Nstart+1)-x(Nstart))
!    Fend = (V(size(V))-V(size(V)-1))/(x(size(x))-x(size(x)-1))
!    a(2) = abs(2./((F2-a(1))/(x(Nstart)-x(Nstart-1))))
!    a(3) = a(1)/Fend
    
!    call mrqmin(x,V,sig,a,mask,covar,alpha,chisq,funcs,alamda)

!    call funcs(x,a,Vfitted,dfit)
    
!    if (present(plot)) then
!        if (plot) then
!            call plt%initialize(grid=.true.,xlabel='$x [nm]$',ylabel='$V (V)$', &
!                            figsize=[20,15], font_size=20, title='Fitting', &
!                            legend=.true.,axis_equal=.false., legend_fontsize=20, &
!                            xtick_labelsize=20,ytick_labelsize=20,axes_labelsize=20)
            
!            call plt%add_plot(x,V,label='$data$', &
!                    linestyle='bo',linewidth=2)
                    
!            call plt%add_plot(x,Vfitted,label='$fitting$', &
!                                linestyle='r-',linewidth=2)
!            call plt%savefig('png/fitting.png', pyfile='python/fitting.py')
!            print *, a
!        endif
!    endif
    
!    F=a(1)
!    R=a(2)
!    gamma=a(3)
    

!    contains

!    subroutine funcs(x,p,yfit,dydp)
!        real(dp), intent(in)    :: x(:), p(:)
!        real(dp), intent(out)   :: yfit(:)
!        real(dp), intent(out)   :: dydp(:,:)
        
!        yfit = (p(1)*p(2)*x*(p(3)-1.d0)+p(1)*x**2) / (p(3)*x+p(2)*(p(3)-1.d0))
!        dydp(:,1) = (p(2) * (p(3)-1.d0)*x + x**2) / &
!                    (p(2) * (p(3)-1.d0) + p(3)*x)
!        dydp(:,2) = (p(1) * (x**2) * (p(3)-1.d0)**2) / &
!                    (p(2) * (p(3)-1.d0) + p(3) * x)**2
!        dydp(:,3) = -(p(1) * x**3) / (p(2) * (p(3)-1.d0) + p(3) * x)**2

!    end subroutine funcs
    
!end subroutine fitpot

subroutine fitpot(x,V,F,R,gamma,plot)
    !Fits V(x) data to F,R,gamma standard potential using L-V
    !minimization module
    use Levenberg_Marquardt, only: nlinfit
    use pyplot_mod, only: pyplot
    
    real(dp), intent(in)    ::x(:),V(:)
    logical, intent(in), optional   :: plot
    real(dp), intent(out)   ::F,R,gamma
    real(dp)                ::p(3),F2,Fend,var
    
    integer                 :: Nstart=3,i
    
    type(pyplot)            :: plt
  
    
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
    
    if (present(plot)) then
        if (plot) then
            call plt%initialize(grid=.true.,xlabel='$x [nm]$',ylabel='$V (V)$', &
                            figsize=[20,15], font_size=20, title='Fitting', &
                            legend=.true.,axis_equal=.false., legend_fontsize=20, &
                            xtick_labelsize=20,ytick_labelsize=20,axes_labelsize=20)
            
            call plt%add_plot(x,V,label='$data$', &
                    linestyle='bo',linewidth=2)
                    
            call plt%add_plot(x,[(fun(x(i),p),i=1,size(x))],label='$fitting$', &
                                linestyle='r-',linewidth=2)
            call plt%savefig('png/fitting.png', pyfile='python/fitting.py')
            
            print *, p
        endif
    endif

    contains

    pure function fun(x,p) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: p(:)
        real(dp)             :: y
        y=(p(1)*p(2)*x*(p(3)-1.d0)+p(1)*x**2) / (p(3)*x+p(2)*(p(3)-1.d0))
    end function fun
end subroutine fitpot

end module interface_helmod