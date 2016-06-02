module heating

implicit none

integer, parameter  :: dp=8, Nr=64
integer, parameter  :: kx=4, ky=4, kz=4, iknot=0  !interpolation parameters
real(dp), parameter :: Jlim=1.d-19, Fmin = 1.8d0 , rho = 16.8d0 !(ohm * nm)




logical, parameter  :: debug = .true.

type, public        :: HeatData

    real(dp)        :: dx=1.d-6, dt = 1.d-5, convergence_criterion = 1.d-15
    integer         :: solvesteps

    real(dp)        :: tempinit(:), tempfinal(:), T_boundary
    real(dp)        :: fse

type, private       :: InterData

    real(dp)                :: F(3) = [1.d0, 1.d0, 1.d0], kT = 5.d-2, W = 4.5d0
    real(dp)                :: Jem, heat !ouput parameters
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
    real(dp)                :: timings(7) !timing variables
    !1: Interpolation set, 2: interpolation evaluate 3: Cur: fitting
    !4: Cur: 1D interpolation set. 
    !5: Cur: J_num_integ with barrier model
    !6: J_num_integ with interpolation
    !7: Cur: GTF
    integer                 :: counts(7)  !counting variables (same as timings)
end type InterData

contains
 
subroutine interp_set(this)
    use bspline, only: db3ink

    type(InterData), intent(inout)      :: this
    real(dp), allocatable               :: x(:), y(:), z(:)
    integer                             :: i, iflag
    real(dp)                            :: t2, t1 !timing
    
    if (debug) call cpu_time(t1)
    this%nx = size(this%phi,1)
    this%ny = size(this%phi,2)
    this%nz = size(this%phi,3)


    allocate(x(this%nx),y(this%ny),z(this%nz),this%bcoef(this%nx,this%ny,this%nz))
    allocate(this%tx(this%nx + kx),this%ty(this%ny + ky),this%tz(this%nz + kz))
    x = [(this%grid_spacing(1)*(i-1), i=1,this%nx)]
    y = [(this%grid_spacing(2)*(i-1), i=1,this%ny)]
    z = [(this%grid_spacing(3)*(i-1), i=1,this%nz)]
    call db3ink(x, this%nx, y, this%ny, z, this%nz, this%phi, kx, ky, kz, iknot, &
                this%tx, this%ty, this%tz, this%bcoef, iflag)
    if (iflag /= 0) then
        stop 'Something went wrong with interpolation set'
    endif
    this%set = .true.
    
    if (debug) then
        call cpu_time(t2)
        this%timings(1) = this%timings(1) + t2 - t1
        this%counts(1) = this%counts(1) + 1
    endif
end subroutine interp_set
 


function get_heat(this) result(TotalHeat)
    use pyplot_mod, only: pyplot

    type(InterData), intent(inout)  :: this
    
    integer             :: i, j, k, Nx, Ny, Nz, Nxy, icount, ptype 
    real(dp)            :: dS
    real(dp), dimension(size(this%phi,3))  :: Pnot, Icur, Ar, TotalHeat , Pjoule
    
    type(pyplot)        :: plt
    integer, parameter  :: font = 35
    
    Nx = size(this%phi,1)
    Ny = size(this%phi,2)
    Nz = size(this%phi,3)
    Nxy=Nx*Ny
    dS = product(this%grid_spacing(1:2))
    
    Icur(Nz) = 0.d0
    do k = Nz-1, 2, -1
        icount=0
        Pnot(k) = 0.d0
        Icur(k) = Icur(k+1) / dS
        Ar(k) = 0.d0
        do j=2,Ny-1
            do i=2,Nx-1
                ptype = classify_point(i,j,k)
                if (ptype == 0) then !surface point
                    this%Nstart = [i,j,k]
                    call J_from_phi(this)
                    Pnot(k) = Pnot(k) + this%heat !add Not heat 
                    Icur(k) = Icur(k) + this%Jem  ! add Current to slice
                    Ar(k) = Ar(k) + dS
                    icount = icount + 1
                else if (ptype == -1) then !interior point
                    Ar(k) = Ar(k) + dS   !increase slice's area
                    icount = icount + 1
                endif
                if (icount > Nx*Ny / 2) exit
            enddo
            if (icount > Nx*Ny / 2) then
                Pnot(k) = 0.d0 
                Icur(k) =0.d0  
                Ar(k) = dS * Nx * Ny
                exit
            endif
        enddo
        Icur(k) = Icur(k) * dS
        Pnot(k) = Pnot(k) * dS
        Pjoule(k) =  rho * this%grid_spacing(3) * Icur(k)**2 / max(Ar(k),1.d-100)
        TotalHeat(k) =  Pjoule(k) + Pnot(k)
        if (icount > 1 .and. icount < Nx*Ny/2) print *, 'icount=',icount,'Icur=', &
            Icur(k), 'Pnot =', Pnot(k), 'Pjoule=', Pjoule(k), 'Ptot=', TotalHeat(k) 
    enddo
    
    call plt%initialize(grid=.true.,xlabel='$z [nm]$',ylabel='$heat [W]$', &
            figsize=[20,15], font_size=font, title='FN-plot test', &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
    call plt%add_plot([(this%grid_spacing(3)*i, i=0,size(TotalHeat)-1)], TotalHeat, &
                label='$current$', linestyle='b-',linewidth=2, yscale = 'log')
                    
    call plt%savefig('png/tipheat.png', pyfile='python/tipheat.py')
    
    contains
    
    function classify_point(i,j,k) result(pointype)
    
        integer, intent(in) :: i,j,k !point indices
        integer             :: pointype !1: exterior, 0: surface, -1: interior
    
        if (this%phi(i,j,k) > 1.d-8) then
            pointype = 1
        else
            if ((this%phi(i+1,j,k)>1.d-8 .or. this%phi(i-1,j,k)>1.d-8 .or. &
                 this%phi(i,j+1,k)>1.d-8 .or. this%phi(i,j-1,k)>1.d-8 .or. &
                 this%phi(i,j,k+1)>1.d-8 .or. this%phi(i,j,k-1)>1.d-8)) then
                pointype = 0
            else
                pointype = -1
            endif
        endif
    end function classify_point

end function get_heat



subroutine J_from_phi(this)
    use bspline, only: db3val
    use std_mat, only: linspace
    use emission, only: EmissionData, cur_dens
    
    type(InterData), intent(inout)  :: this
    
    type(EmissionData),save :: that !emission data structure
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    !the above are spline-related parameters
    
    integer                 :: i, istart,jstart,kstart   
    real(dp)                :: direc(3), Fnorm, rmax, Vmax, xmax, ymax, zmax
                                        !rmax: maximum distance for rline
    logical                 :: badcontition
    real(dp)                :: t1, t2

    badcontition = .false.
    if (debug) call cpu_time(t1)
    
    istart = this%Nstart(1)
    jstart = this%Nstart(2)
    kstart = this%Nstart(3)
    
    !find direction of the line (same as Efield direction)
    this%F = ([this%phi(istart+1,jstart,kstart), &
                 this%phi(istart,jstart+1,kstart), &
                 this%phi(istart,jstart,kstart+1)] - &
                 [this%phi(istart-1,jstart,kstart), &
                 this%phi(istart,jstart-1,kstart), &
                 this%phi(istart,jstart,kstart-1)])/this%grid_spacing
    Fnorm = norm2(this%F)
    if (Fnorm> Fmin) then
        direc = this%F / Fnorm
        rmax = 2.d0 * this%W / Fnorm
        do i = 1,10
            xmax = this%grid_spacing(1)*(istart-1) + rmax * direc(1)
            ymax = this%grid_spacing(2)*(jstart-1) + rmax * direc(2)
            zmax = this%grid_spacing(3)*(kstart-1) + rmax * direc(3)
            call db3val(xmax, ymax, zmax, idx, idy, idz, &
                    this%tx, this%ty, this%tz, this%nx, this%ny, this%nz, &
                    kx, ky, kz, this%bcoef, Vmax,  &
                    iflag, inbvx, inbvy, inbvz, iloy, iloz)
            if (iflag /=0) then
                rmax = rmax / 1.3d0
                badcontition = .true.
                exit
            endif
            
            if (Vmax > this%W) exit
            rmax = 1.3d0 * rmax 
        enddo
        if (i==11) badcontition = .true.
        
        this%rline = linspace(0.d0, rmax, Nr)
        !set the line of interpolation and interpolate
        this%xline = this%grid_spacing(1)*(istart-1) + this%rline * direc(1)
        this%yline = this%grid_spacing(2)*(jstart-1) + this%rline * direc(2)
        this%zline = this%grid_spacing(3)*(kstart-1) + this%rline * direc(3)
       
        do i=1,Nr !interpolate for all points of rline
            call db3val(this%xline(i), this%yline(i), this%zline(i), idx, idy, idz, &
                        this%tx, this%ty, this%tz, this%nx, this%ny, this%nz, &
                        kx, ky, kz, this%bcoef, this%Vline(i),  &
                        iflag, inbvx, inbvy, inbvz, iloy, iloz)
        enddo
        that%xr = this%rline
        that%Vr = this%Vline
        that%mode = 1 !in general case use interpolation
        if (badcontition) that%mode = -3 
    else
        that%F = norm2(this%F)
        that%R = 1.d4
        that%gamma = 1.d0
        that%mode = 0
    endif
    
    if (debug) then
        call cpu_time(t2)
        this%timings(2) = this%timings(2) + t2 - t1
        this%counts(2) = this%counts(2) + 1
    endif
    
    that%W = this%W
    that%kT = this%kT
    that%full = .false.
    
    call cur_dens(that)
    
    if (that%Jem > Jlim) then
        that%full = .true.
        that%mode = -2 !try fitting first. If not successful, do interpolation
        call cur_dens(that)
    endif
    this%Jem = that%Jem
    this%heat = that%heat
    do i =1,5 !copy timings and counts from that
        this%timings(i+2) = that%timings(i)
        this%counts(i+2) = that%counts(i)
    enddo

end subroutine J_from_phi


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

end module new_interface
