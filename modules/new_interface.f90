module interface_helmod

implicit none

integer, parameter  :: dp=8, Nr=32 
integer, parameter  :: kx=2, ky=2, kz=2, iknot=0  !interpolation parameters
real(dp), parameter :: Jlim=1.d-22, Fmin = 0.5d0, rmax = 3.d0

type, public       :: InterData

    real(dp)                :: F(3), kT, W, Jem, heat
    !Parameters defining emission at each point
    real(dp)                :: grid_spacing(3)
    !External potential and grid spacing parameters
    integer                 :: Nstart(3), nx, ny, nz
    !indices of starting surface surf_points
    real(dp)                :: rline(Nr), Vline(Nr), xline(Nr), yline(Nr), zline(Nr)
    !line data values -- potential
    
    real(dp), allocatable   :: phi(:,:,:), bcoef(:,:,:), tx(:), ty(:), tz(:)
    !External potential and its corresponding spline parameters
    logical                 :: set = .false.
    !true if interpolation has been set (db3ink called)
    real(dp)                :: TimeCur=0.d0, TimeInSet=0.d0, TimeInt=0.d0 
    !timing variables
end type InterData
    

contains
 
subroutine interp_set(this)
    use bspline, only: db3ink

    type(InterData), intent(inout)      :: this
    real(dp), allocatable               :: x(:), y(:), z(:)
    
    this%nx = size(this%phi,1)
    this%ny = size(this%phi,2)
    this%nz = size(this%phi,3)


    allocate(x(this%nx),y(this%ny),z(nthis%z),this%bcoef(this%nx,this%ny,this%nz))
    allocate(this%tx(this%nx + kx),this%ty(this%ny + ky),this%tz(this%nz + kz))
    x = [(this%grid_spacing(1)*(i-1), i=1,nx)]
    y = [(this%grid_spacing(2)*(i-1), i=1,ny)]
    z = [(this%grid_spacing(3)*(i-1), i=1,nz)]
    call db3ink(x, this%nx, y, this%ny, z, this%nz, this%phi, kx, ky, kz, iknot, &
                this%tx, this%ty, this%tz, this%bcoef, iflag)
    if (iflag /= 0) then
        stop 'Something went wrong with interpolation set'
    endif
    
end subroutine interp_set
 
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

subroutine J_from_phi(this)
    use bspline, only: db3val
    use std_mat, only: linspace
    use emission, only: EmissionData, cur_dens
    
    type(InterData), intent(inout)  :: this
    
    type(EmissionData)      :: that !emission data structure
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    !the above are spline-related parameters
    
    integer                 :: i, j, istart,jstart,kstart
    real(dp), dimension(3)  :: direc, Efstart
    
    logical                 :: plot = .false.
    
    real(dp)                :: t1, t2

    
    call cpu_time(t1)

    this%rline = linspace(0.d0, rmax, Nr)

    istart = this%Nstart(1,j)
    jstart = this%Nstart(2,j)
    kstart = this%Nstart(3,j)
    
        !find direction of the line (same as Efield direction)
        ! always 
    this%F = ([this%phi(istart+1,jstart,kstart), &
                 this%phi(istart,jstart+1,kstart), &
                 this%phi(istart,jstart,kstart+1)] - &
                 [this%phi(istart-1,jstart,kstart), &
                 this%phi(istart,jstart-1,kstart), &
                 this%phi(istart,jstart,kstart-1)])/this%grid_spacing
                 
    if (norm2(this%F)> Fmin) then
        direc = this%F / norm2(this%F)
        !set the line of interpolation and interpolate
        this%xline = x(istart) + this%rline * direc(1)
        this%yline = y(jstart) + this%rline * direc(2)
        this%zline = z(kstart) + this%rline * direc(3)
                
        do i=1,Nr !interpolate for all points of rline
            call db3val(this%xline(i), this%yline(i), this%zline(i), idx, idy, idz, &
                        this%tx, this%ty, this%tz, this%nx, this%ny, this%nz, &
                        kx, ky, kz, this%bcoef, this%Vline(i),  &
                        iflag, inbvx, inbvy, inbvz, iloy, iloz)
        enddo
        call cpu_time(t2)
        this%TimeInt = this%TimeInt + t2 - t1
        
        that%xr = this%rline
        that%Vr = this%Vline
        that%mode = 1
    else
        that%F = norm2(Efstart)
        that%R = 1.d4
        that%gamma = 1.d0
        that%mode = 0
    endif
    
    that%W = this%W
    that%kT = this%kT
    
    call cpu_time(t1)
    call cur_dens(that)
    if (that%Jem > Jlim) then
        that%full = .true.
        call cur_dens(that)
    endif
    call cpu_time(t2)
    this%TimeCur = this%TimeCur + t2 - t1

end subroutine J_from_phi


