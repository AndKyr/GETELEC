module heating

implicit none

integer, parameter  :: dp=8, Nr=64
real(dp), parameter :: kBoltz=8.6173324d-5
integer, parameter  :: kx=4, ky=4, kz=4, iknot=0  !interpolation parameters
real(dp), parameter :: Jlim=1.d-19, Fmin = 1.8d0 , rho = 16.8d0 !(ohm * nm)
real(dp), parameter :: dt = 1.d-5, convergence_criterion = 1.d-15, workfunc = 4.5d0

real(dp)            :: fse = 3.d0



logical, parameter  :: debug = .true.

!type, public        :: HeatData

!    integer         :: solvesteps  = 100000, tipheight = 0

!    real(dp)        :: tempinit(:), tempfinal(:), heat(:), T_boundary
   

!end type HeatData

type, public        :: Potential !data defining the potential at all gridpoints

    real(dp), allocatable   :: phi(:,:,:) !potential at all gridpoints
    real(dp)                :: grid_spacing(3)
    integer                 :: nx, ny, nz !shape(phi)
    
    real(dp),allocatable    :: bcoef(:,:,:), tx(:), ty(:), tz(:) 
    !interpolation data. Set after db3ink is called
    logical                 :: set = .false.
    !true if interpolation has been set (db3ink called)
    real(dp)                :: timing
end type Potential

type, public        :: PointEmission
   
    real(dp)                :: F(3) = [1.d0, 1.d0, 1.d0], kT = 5.d-2, W = 4.5d0
    real(dp)                :: Jem, heatNot !ouput parameters
    !Parameters defining emission at each point

    integer                 :: Nstart(3)
    !indices of starting surface surf_points
    real(dp)                :: rline(Nr), Vline(Nr), xline(Nr), yline(Nr), zline(Nr)
    !line data values -- potential
    
    real(dp)                :: timings(7) !timing variables
    !1: Interpolation set, 2: interpolation evaluate 3: Cur: fitting
    !4: Cur: 1D interpolation set. 
    !5: Cur: J_num_integ with barrier model
    !6: J_num_integ with interpolation
    !7: Cur: GTF
    integer                 :: counts(7)  !counting variables (same as timings)
end type PointEmission

contains

function get_heat(temp,poten) result(TotalHeat)
    
    type(Potential), intent(in)     :: poten !inputed potential
    real(dp), intent(in)            :: temp(:) !temperature distribution along tip
    
    type(PointEmission)             :: point
               
    integer                         :: i, j, k, icount, ptype !type of point 
    real(dp)                        :: dS, rho 
    !area corresponding at each gridpoint, resistivity
    real(dp), dimension(size(temp)) :: Pnot, Icur, Ar, TotalHeat , Pjoule
    !Not heat (W), Total Current through slice(Amps), Slice area (nm**2),
    !total heat density at slice(W / m^3), Joule heat at lattice(W)                                             
    dS = product(poten%grid_spacing(1:2))
    point%W = workfunc
    Icur = 0.d0
    do k = poten%nz-1, 2, -1
        icount=0
        Pnot(k) = 0.d0 !Not counts from zero
        Icur(k) = Icur(k+1) / dS !upwards current is accumulative
        Ar(k) = 1.d-20
        point%kT = kBoltz * temp(k)
        rho = resistivity(temp(k)) * fse !resistivity of the slice 
        do j=2,poten%ny-1
            do i=2,poten%nx-1
                ptype = classify_point(i,j,k)
                if (ptype == 0) then !surface point
                    point%Nstart = [i,j,k]
                    call emit_atpoint(poten,point)
                    Pnot(k) = Pnot(k) + point%heatNot !add Not heat 
                    Icur(k) = Icur(k) + point%Jem  ! add Current to slice
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
        Icur(k) = Icur(k) * dS
        Pnot(k) = Pnot(k) * dS
        Pjoule(k) =  rho * poten%grid_spacing(3) * Icur(k)**2 / max(Ar(k), 1.d-100)
        TotalHeat(k) =  1.d27*(Pjoule(k) + Pnot(k)) / (Ar(k) * poten%grid_spacing(3))
        !Total heat density in SI units
        if (icount > 1 .and. icount < poten%Nx * poten%Ny/2) &
            print *, 'icount=',icount,'Icur=', Icur(k), 'Pnot =', Pnot(k), &
                    'Pjoule=', Pjoule(k), 'Ptot=', TotalHeat(k) 
    enddo
    
    contains
    
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

end function get_heat

subroutine emit_atpoint(poten,point)
    use bspline, only: db3val
    use std_mat, only: linspace
    use emission, only: EmissionData, cur_dens
    
    type(Potential), intent(in)         :: poten !potential data struct
    type(PointEmission), intent(inout)  :: point !point data struct
    
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
    
    if (.not. poten%set) print *, 'Error: interpolation not set'
    
    if (debug) call cpu_time(t1)
    istart = point%Nstart(1)
    jstart = point%Nstart(2)
    kstart = point%Nstart(3)
    
    !find direction of the line (same as Efield direction)
    point%F = ([poten%phi(istart+1,jstart,kstart), &
                 poten%phi(istart,jstart+1,kstart), &
                 poten%phi(istart,jstart,kstart+1)] - &
                 [poten%phi(istart-1,jstart,kstart), &
                 poten%phi(istart,jstart-1,kstart), &
                 poten%phi(istart,jstart,kstart-1)])/poten%grid_spacing
    Fnorm = norm2(point%F)
    if (Fnorm> Fmin) then
        direc = point%F / Fnorm
        rmax = 2.d0 * point%W / Fnorm
        do i = 1,10
            xmax = poten%grid_spacing(1)*(istart-1) + rmax * direc(1)
            ymax = poten%grid_spacing(2)*(jstart-1) + rmax * direc(2)
            zmax = poten%grid_spacing(3)*(kstart-1) + rmax * direc(3)
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
        point%xline = poten%grid_spacing(1)*(istart-1) + point%rline * direc(1)
        point%yline = poten%grid_spacing(2)*(jstart-1) + point%rline * direc(2)
        point%zline = poten%grid_spacing(3)*(kstart-1) + point%rline * direc(3)
       
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
    
    if (that%Jem > Jlim) then
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

end subroutine emit_atpoint


subroutine dsolve(T, heat, Tb, dx, dt)
    implicit none
    double precision, intent(inout) :: T(:) !< Temperature distribution of the tip. In: previous estimate, out: current estimate
    double precision, intent(in) :: heat(:) !< Deposited heat distribution as a function of height
    double precision, intent(in) :: Tb !< Boundary temperature at the bottom of the tip
    double precision, intent(in) :: dx !< Grid size (m)
    double precision, intent(in) :: dt !< Timestep (s)

    double precision :: Told(size(T)) ! Old temperature distribution
    double precision, parameter :: Cv = 3.45d6 ! J/(m^3 K)
    double precision :: K ! Thermal conductivity, W/(m K)
    double precision :: K1 
    double precision, parameter :: K2 = 1d0/Cv ! rho0/(Cv*T0)
    double precision :: K3
    double precision :: a(size(T)), b(size(T)), c(size(T)), d(size(T))
    double precision :: cd
    integer :: i
    double precision, parameter :: L = 2.44e-8 ! Wiedemann Frantz
    logical, save :: firsttime = .true.
    integer :: height

    integer, parameter :: method = 2

    height = size(T)
    if(firsttime) then
        print *, "fse", fse, "method", method
        print *, "J(height)", J(height)
        firsttime = .false.
    end if

    rho = 17e-9*fse ! Electrical resistivity
    K = 400d0/fse ! Thermal conductivity

    K1 = K/Cv

    Told = T
    T(1) = Tb

    if(method==1) then
      ! Schuster et. al resistivity
        do i = 2, height-1
            rho = resistivity(T(i))*1d-8*fse
            K = L*T(i)/rho
            K1 = K/Cv

            a(i) = -K1*dt/(2*dx**2)
            b(i) = 1+K1*dt/dx**2
            c(i) = -K1*dt/(2*dx**2)
            ! Here we cheat a bit. Use old value of T as temperature when it should be avg. of new and old,
        ! since the resistivity is non-linear (Can't solve for the new temperature
        ! otherwise)
        d(i) = K1*dt/(2*dx**2)*(Told(i-1)-2*Told(i)+Told(i+1))+K2*dt*heat
        end do

        a(1) = 0
      ! Boundry conditions
      ! T(1) = Tb =>
        b(1) = 1d0
        c(1) = 0d0
        d(1) = Tb

      ! dT/dx|height = 0 => T(height) = T(height-1) =>
        a(height) = 1d0
        b(height) = -1d0
        d(height) = 0d0

      ! Solve tridiagonal matrix with above coefficients
        call tridiag(a, b, c, d, T)
    else if (method==2) then
        K3 = J(height)**2*K2
      ! Backward Euler with simple linear resistivity. Useful when e.g. debugging
        a = -K1*dt/dx**2
        b = 1+2*K1*dt/dx**2
        c = -K1*dt/dx**2
        d = T+K2*heat*dt

      ! Boundry conditions
        b(1) = 1d0
        c(1) = 0d0
        d(1) = Tb

        a(height) = 1d0
        b(height) = -1d0
        d(height) = 0d0

      ! Solve tridiagonal matrix with above coefficients
        call tridiag(a, b, c, d, T)
    else if (method==3) then
      ! Very simple Forward Euler. If using this you still get incorrect results,
      ! clearly the problem must not be in the solver, but somewhere else.
        T(1) = Tb
        Told = T
        do i = 2, height-1
            cd = (Told(i-1)-2*Told(i)+Told(i+1))/dx**2
            T(i) = (K1*cd+K2*heat(i))*dt + Told(i)
        end do
        i = height

        cd = (Told(i-1)-2*Told(i)+Told(i-1))/dx**2
        T(i) = (K1*cd+K2*heat(i))*dt + Told(i)
    else
        stop "Unknown solving method"
    end if

    T(1) = Tb
end subroutine

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

end module heating
