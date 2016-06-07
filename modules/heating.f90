module heating

implicit none

integer, parameter  :: dp = 8, Nr = 64 !No of interpolation points in emission line
real(dp), parameter :: kBoltz = 8.6173324d-5 !Boltzmann constant in eV/K
integer, parameter  :: kx = 4, ky = 4, kz = 4, iknot = 0  !interpolation parameters
real(dp), parameter :: Jlimratio = 1.d-4, Flimratio = 0.2d0 
!minimum current and field ratio to max for full calculation
real(dp), parameter :: convergence_criterion = 1.d-15, workfunc = 4.5d0

real(dp)            :: fse = 50.d0 !finite size effect for conductivity

logical, parameter  :: debug = .true.

type, public        :: HeatData

    real(dp), allocatable   :: tempinit(:), tempfinal(:), hpower(:), J_avg(:)
    !initial temperature, final temperature, deposited heating power (W/nm^3)
    real(dp)                :: Tbound, maxtime, dt = 1.d-3, dx !dx in nm
    !boundary bulk temperature, maximum time evolution, timestep, spacestep, all SI
    integer                 :: tipbounds(2), maxsteps = 1e6, stepsdone
    ! tipbounds in grid, maximum time steps
    real(dp)                :: temperror
end type HeatData

type, public        :: Potential !data defining the potential at all gridpoints

    real(dp), allocatable   :: phi(:,:,:) !potential at all gridpoints
    real(dp)                :: grid_spacing(3)  !in nm
    integer                 :: nx, ny, nz !shape(phi)
    
    real(dp),allocatable    :: bcoef(:,:,:), tx(:), ty(:), tz(:) 
    !interpolation data. Set after db3ink is called
    logical                 :: set = .false.
    !true if interpolation has been set (db3ink called)
    real(dp)                :: timing
end type Potential

type, public        :: PointEmission
   
    real(dp)                :: F(3) = [1.d0, 1.d0, 1.d0], kT = 5.d-2, W = 4.5d0
    !F in V/nm, kT, W in eV
    real(dp)                :: Jem, heatNot !ouput parameters
    !Parameters defining emission at each point. Jem in A/nm^2
    real(dp)                :: Jmax = 0.d0, Fmax = 0.d0 
    !maximum current,field encountered up to now

    integer                 :: Nstart(3)
    !indices of starting surface surf_points
    real(dp)                :: rline(Nr), Vline(Nr), xline(Nr), yline(Nr), zline(Nr)
    !line data values -- potential . All in nm
    
    real(dp)                :: timings(7) !timing variables
    !1: Interpolation set, 2: interpolation evaluate 3: Cur: fitting
    !4: Cur: 1D interpolation set. 
    !5: Cur: J_num_integ with barrier model
    !6: J_num_integ with interpolation
    !7: Cur: GTF
    integer                 :: counts(7)  !counting variables (same as timings)
end type PointEmission

contains

subroutine get_heat(heat,poten)
! finds surfacepoints from poten, calculates current and Nottingham at each one
! and calculates deposited heat density at each grid z slice    
    type(Potential), intent(in)     :: poten !inputed potential
    type(HeatData), intent(inout)   :: heat !heat data inout
        
    type(PointEmission)             :: point
               
    integer                         :: i, j, k, icount, ptype !type of point 
    real(dp)                        :: dS, rho 
    !area corresponding at each gridpoint, resistivity
    real(dp), dimension(size(heat%tempinit)) :: Pnot, Icur, Ar, Pjoule
    !Not heat (W), Total Current through slice(Amps), Slice area (nm**2),
    !total heat density at slice(W / m^3), Joule heat at lattice(W)
    logical                         :: intip
                                                 
    dS = product(poten%grid_spacing(1:2))
    point%W = workfunc
    Icur = 0.d0
    intip = .false.
    heat%dx = poten%grid_spacing(3)
    do k = poten%nz-1, 2, -1
        icount=0
        Pnot(k) = 0.d0 !Not counts from zero
        Icur(k) = Icur(k+1) / dS !upwards current is accumulative
        Ar(k) = 1.d-20
        point%kT = kBoltz * heat%tempinit(k)
        rho = resistivity(heat%tempinit(k)) * fse !resistivity of the slice 
        do j=2,poten%ny-1
            do i=2,poten%nx-1
                ptype = classify_point(i,j,k)
                if (ptype == 0) then !surface point
                    point%Nstart = [i,j,k]
                    call emit_atpoint(poten,point)
                    Pnot(k) = Pnot(k) + point%heatNot !add Not heat 
                    Icur(k) = Icur(k) + point%Jem  ! add Current density to slice
                    Ar(k) = Ar(k) + dS
                    icount = icount + 1
                else if (ptype == -1) then !interior point
                    Ar(k) = Ar(k) + dS   !increase slice's area
                    icount = icount + 1
                endif
                if (icount > poten%Nx * poten%Ny / 2) exit !not at tip
            enddo
            if (icount > poten%Nx * poten%Ny / 2) then  !not at tip
                Pnot(k) = 0.d0 
                Icur(k) =0.d0  
                Ar(k) = dS * poten%Nx * poten%Ny
                exit
            endif
        enddo
        Icur(k) = Icur(k) * dS !multiply by elementary area to get current
        Pnot(k) = Pnot(k) * dS !multiply by elementary area to get Watts
        Pjoule(k) =  rho * heat%dx * (Icur(k)**2) / max(Ar(k), 1.d-100)
        heat%hpower(k) =  (Pjoule(k) + Pnot(k)) / (Ar(k) * heat%dx)
        !Total heat density in SI units
        if ((.not. intip) .and. (Icur(k) >0.d0)) then
            heat%tipbounds(2) = k
            intip = .true.
        endif
        if (intip .and. Icur(k) <= 0.d0) then
            heat%tipbounds(1) = k + 1
            intip = .false.
        endif
        if (debug) then
            if (intip) print *, 'icount=',icount,'Icur=', Icur(k), 'Pnot =', &
                Pnot(k), 'Pjoule=', Pjoule(k), 'Ptot=', heat%hpower(k) 
        endif
    enddo
    heat%J_avg = Icur / Ar
    
    if (debug) print *, 'Timings', point%timings
    
    contains
    
    function classify_point(i,j,k) result(pointype)
    
        integer, intent(in) :: i,j,k !point indices
        integer             :: pointype !1: exterior, 0: surface, -1: interior
    
        if (poten%phi(i,j,k) > 1.d-8) then
            pointype = 1
        else
            if ((poten%phi(i+1,j,k)>1.d-8 .or. poten%phi(i-1,j,k)>1.d-8 .or. &
                 poten%phi(i,j+1,k)>1.d-8 .or. poten%phi(i,j-1,k)>1.d-8 .or. &
                 poten%phi(i,j,k+1)>1.d-8 .or. poten%phi(i,j,k-1)>1.d-8)) then
                pointype = 0
            else
                pointype = -1
            endif
        endif
    end function classify_point

end subroutine get_heat

subroutine emit_atpoint(poten,point)
! calculates emission current and nottingham heating in a grid point
! carefull: every distance in this function is in nm 
    use bspline, only: db3val
    use std_mat, only: linspace
    use GeTElEC, only: EmissionData, cur_dens
    
    type(Potential), intent(in)         :: poten !potential data struct
    type(PointEmission), intent(inout)  :: point !point data struct
    
    type(EmissionData),save :: that !emission data structure
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    !the above are spline-related parameters
    
    integer                 :: i, istart,jstart,kstart   
    real(dp)                :: direc(3), Fnorm, rmax, Vmax, xmax, ymax, zmax 
    !direc: direction of field, rmax: maximum distance for rline. All in nm
    real(dp)                :: dx, dy, dz !grid spacing in nm
    logical                 :: badcontition
    real(dp)                :: t1, t2

    badcontition = .false.
    
    if (.not. poten%set) print *, 'Error: interpolation not set'
    
    if (debug) call cpu_time(t1)
    istart = point%Nstart(1)
    jstart = point%Nstart(2)
    kstart = point%Nstart(3)
    dx = poten%grid_spacing(1)
    dy = poten%grid_spacing(2)
    dz = poten%grid_spacing(3)
    
    !find direction of the line (same as Efield direction)
    point%F = ([poten%phi(istart+1,jstart,kstart), &
                 poten%phi(istart,jstart+1,kstart), &
                 poten%phi(istart,jstart,kstart+1)] - &
                 [poten%phi(istart-1,jstart,kstart), &
                 poten%phi(istart,jstart-1,kstart), &
                 poten%phi(istart,jstart,kstart-1)]) / [dx, dy, dz]
    Fnorm = norm2(point%F)
    if (Fnorm > point%Fmax) point%Fmax = Fnorm
    if (Fnorm >= point%Fmax * Flimratio) then
        direc = point%F / Fnorm
        rmax = 2.d0 * point%W / Fnorm
        do i = 1,10
            xmax = dx * (istart-1) + rmax * direc(1)
            ymax = dy * (jstart-1) + rmax * direc(2)
            zmax = dz * (kstart-1) + rmax * direc(3)
            call db3val(xmax, ymax, zmax, idx, idy, idz, &
                    poten%tx, poten%ty, poten%tz, poten%nx, poten%ny, poten%nz, &
                    kx, ky, kz, poten%bcoef, Vmax,  &
                    iflag, inbvx, inbvy, inbvz, iloy, iloz)
            if (iflag /=0) then
                rmax = rmax / 1.3d0
                badcontition = .true.
                exit
            endif
            
            if (Vmax > point%W) exit
            rmax = 1.3d0 * rmax 
        enddo
        if (i==11) badcontition = .true.
        
        point%rline = linspace(0.d0, rmax, Nr)
        !set the line of interpolation and interpolate
        point%xline = dx * (istart-1) + point%rline * direc(1)
        point%yline = dy * (jstart-1) + point%rline * direc(2)
        point%zline = dz * (kstart-1) + point%rline * direc(3)
       
        do i=1,Nr !interpolate for all points of rline
            call db3val(point%xline(i), point%yline(i), point%zline(i), idx, idy, idz, &
                        poten%tx, poten%ty, poten%tz, poten%nx, poten%ny, poten%nz, &
                        kx, ky, kz, poten%bcoef, point%Vline(i),  &
                        iflag, inbvx, inbvy, inbvz, iloy, iloz)
        enddo
        that%xr = point%rline
        that%Vr = point%Vline
        that%mode = 1 !in general case use interpolation
        if (badcontition) that%mode = -3 
    else
        that%F = norm2(point%F)
        that%R = 1.d4
        that%gamma = 1.d0
        that%mode = 0
    endif
    
    if (debug) then
        call cpu_time(t2)
        point%timings(2) = point%timings(2) + t2 - t1
        point%counts(2) = point%counts(2) + 1
    endif
    
    that%W = point%W
    that%kT = point%kT
    that%full = .false.
    
    call cur_dens(that)
    
    if (that%Jem > point%Jmax * Jlimratio) then
        that%full = .true.
        that%mode = -2 !try fitting first. If not successful, do interpolation
        call cur_dens(that)
    endif
    point%Jem = that%Jem
    point%heatNot = that%heat
    do i =1,5 !copy timings and counts from that
        point%timings(i+2) = that%timings(i)
        point%counts(i+2) = that%counts(i)
    enddo
    point%timings(1) = poten%timing
    point%counts(1) = 1
    if (point%Jem > point%Jmax) point%Jmax = point%Jem

end subroutine emit_atpoint


subroutine heateq(heat)
    implicit none

    type(HeatData), intent(inout)   :: heat

    integer, parameter      :: method = 1
    real(dp), parameter     :: Cv = 3.45d-6, L = 2.44e-8 
    ! Cv in(W*fs/(nm^3 K)), L=Wiedemann-Frantz constant in W*Ohm/K^2
    real(dp)                :: K, K1, K2, K3, cd, rho 
    !K: Thermal conductivity, W/(nm K)
    ! K2 = rho0/(Cv*T0)
    real(dp), dimension(heat%tipbounds(2)-heat%tipbounds(1) + 1) :: &
                                                            a, b, c, d, Told, Tnew
    !tridiagonal vectors and temporary temperature working space
    integer :: i, j, solvesteps, Nx
    logical, save :: firsttime = .true.
    
    Nx = heat%tipbounds(2)-heat%tipbounds(1) + 1

    rho = 17.d0*fse ! Electrical resistivity
    K = 3.99e-7/fse ! Thermal conductivity
    solvesteps = min(heat%maxsteps, floor(heat%maxtime / heat%dt))

    K1 = K / Cv
    Tnew = heat%tempinit(heat%tipbounds(1) : heat%tipbounds(2))
    Tnew(1) = heat%Tbound
    
    print *, 'dx = ', heat%dx, 'dt = ', heat%dt
   

    do j = 1, solvesteps
        Told = Tnew
        do i = 2, Nx - 1
            rho = resistivity(Told(i)) * fse
            K = L * Told(i) / rho
            K1 = K / Cv
            a(i) = -K1 * heat%dt / (2.d0 * heat%dx**2)
            b(i) = 1.d0 + K1 * heat%dt / heat%dx**2
            c(i) = -K1 * heat%dt / (2 * heat%dx**2)
            ! Here we cheat a bit. Use old value of T as temperature when it 
            ! should be avg. of new and old, since the resistivity is non-linear.
            ! (Can't solve for the new temperature otherwise)
            d(i) = Told(i) + heat%dt * heat%hpower(i+heat%tipbounds(1)-1) / Cv
        enddo
            
        ! Boundry conditions Tbottom = Tbound
        a(1) = 0.d0
        b(1) = 1.d0
        c(1) = 0.d0
        d(1) = heat%Tbound

        ! dT/dx|height = 0 => T(height) = T(height-1) =>
        a(Nx) = 1.d0
        b(Nx) = -1.d0
        d(Nx) = 0.d0

        ! Solve tridiagonal matrix with above coefficients
        call tridiag(a, b, c, d, Tnew)
        Tnew(1) = heat%Tbound
        heat%temperror = norm2((Tnew - Told) / Told)
        if (heat%temperror < convergence_criterion * heat%dt) exit
    enddo
    heat%tempfinal(heat%tipbounds(1) : heat%tipbounds(2)) = Tnew
    heat%tempfinal(1 : heat%tipbounds(1)) = heat%Tbound
    heat%tempfinal(heat%tipbounds(2):) = -1.d0 
    
    if (debug) then
        print *, 'done ', j-1, 'steps and time = ', j * heat%dt
        print *, 'error = ', heat%temperror
    endif
  
    contains
    
    subroutine tridiag(a, b, c, r, u)
        
        real(dp), intent(in)    :: a(:), b(:), c(:), r(:)
        real(dp), intent(out)   :: u(:)
        integer                 :: n ! Number of rows

        integer                 :: i
        real(dp)                :: b2, c2(size(a))

        n = size(a)

        b2 = b(1)
        u(1) = r(1)/b2
        do i = 2, n
            c2(i) = c(i-1)/b2
            b2 = b(i)-a(i)*c2(i)
            u(i) = (r(i)-a(i)*u(i-1))/b2
        end do

        do i = n-1, 1, -1
            u(i) = u(i)-c2(i+1)*u(i+1)
        end do
    end subroutine tridiag
    
end subroutine heateq

subroutine interp_set(this)
    use bspline, only: db3ink

    type(Potential), intent(inout)      :: this
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
        this%timing = this%timing + t2 - t1
    endif
end subroutine interp_set

function resistivity(T) result(rho)
    real(dp), intent(in)    :: T !temperature
    real(dp)                :: rho
    integer :: Ti
    include "rho_table.h"

    Ti = nint(T)/10

    if (Ti<5 .or. Ti+1>120) then
        print *, "WARNING: Temperature out of range in resisitivity calculation."
        print *, "Temperature capped at 1200K. Your system may be melting!", T

        if(Ti<5) Ti = 5
        if(Ti+1>120) Ti = 120
    end if
    ! Resistivity is an exponential function, so use log to make it linear, 
    !do normal interpolation, and then make it exponential again
    rho = exp( (log(rho_table(Ti+1)) - log(rho_table(Ti))) / (log(Ti*10d0+10)&
                - log(Ti*10d0)) * (log(T) - log(Ti*10d0)) + log(rho_table(Ti)) )

    rho = rho*10.d0 ! Convert to units Ohm*nm
end function resistivity

end module heating