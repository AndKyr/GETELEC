module heating

!module that connects GETELEC to HELMOD. Being input the potential on the grid, it 
!gives the heating at each emitting point on a column-like tip. It also solves the
!heat equation in 1-D.

implicit none

!parameters that are set probably permanently and are not changed by user:

integer, parameter  :: dp = 8, Nr = 64 !No of interpolation points in emission line
real(dp), parameter :: kBoltz = 8.6173324d-5 !Boltzmann constant in eV/K
integer, parameter  :: kx = 4, ky = 4, kz = 4, iknot = 0  !interpolation parameters
real(dp), parameter :: Jlimratio = 1.d-4, Flimratio = 0.2d0 
!minimum current and field ratio to max for full calculation
real(dp), parameter :: convergence_criterion = 1.d-15
!convergence criterion for heat equation in case the steady state is asked
logical, parameter  :: debug = .false., timings = .false., printheat = .true., &
                        savedata = .true. !choose if heating data will be saved
                        
integer, parameter  :: Nsavedata = 50 !save data every Nsavedata steps



!parameters changed by the user that and will be read from file on execution:

real(dp), save      :: workfunc = 4.5d0, fse = 40.d0, tau = 500.d0
!work function, finite size effect multiplyier, brendsen scaling velocity relax time
integer, save       :: compmode = 1, Nsteps
!mode of calculation for comparison purposes
!1: full calculation with full implementation of getelec
!2: forcing "blunt" calculation with GTF (takes into account Nottingham)
!3: Not taking into account the Nottingham effect (only Joule heating)
!4: both (2) and (3)



!Definition of module types
type, public        :: HeatData

    real(dp), allocatable       :: tempinit(:), tempfinal(:), hpower(:), J_avg(:)
    !initial temperature, final temperature, deposited heating power (W/nm^3)
    real(dp)                    :: Tbound, maxtime, dt = 1.d-3, dx 
    !boundary bulk temperature (K), maximum time evolution (fs),
    !timestep (fs), spacestep (nm)
    integer                     :: tipbounds(2), maxsteps = 1e6, stepsdone
    ! tipbounds in grid, maximum time steps
    real(dp)                    :: temperror
end type HeatData

type, public        :: Potential !data defining the potential at all gridpoints

    real(dp), allocatable   :: phi(:,:,:) !potential at all gridpoints
    real(dp)                :: grid_spacing(3)  !in nm
    integer                 :: nx = 0, ny = 0, nz = 0 !shape(phi)
    
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
    
    integer                 :: ierr !error indicator
    !0 : everything went fine
    !10+ getelec ierr : geterlec error in the first not full calculation
    !20+ getelec ierr : getelec error in the second full calculation
    
end type PointEmission


contains


!Module subroutines...
subroutine get_heat(heat,poten)
! finds surfacepoints from poten, calculates current and Nottingham at each one
! and calculates deposited heat density at each grid z slice    
    type(Potential), intent(in)     :: poten !inputed potential
    type(HeatData), intent(inout)   :: heat !heat data inout
        
    type(PointEmission)             :: point
    
    integer, parameter              :: fid = 897465
               
    integer                         :: i, j, k, icount, ptype, Nerrors, Ncorrect 
    real(dp)                        :: dS, rho 
    !area corresponding at each gridpoint, resistivity
    real(dp), dimension(size(heat%tempinit)) :: Pnot, Icur, Ar, Pjoule
    !Not heat (W), Total Current through slice(Amps), Slice area (nm**2),
    !total heat density at slice(W / m^3), Joule heat at lattice(W)
    logical                         :: intip
    integer, save                   :: Ncalls = 0

                                                     
    dS = product(poten%grid_spacing(1:2))
    point%W = workfunc
    Icur = 0.d0
    intip = .false.
    heat%dx = poten%grid_spacing(3)
    point%Fmax = 0.d0
    point%Jmax = 0.d0
    Nerrors = 0 !initialize to zero the number of encountered errors
    Ncorrect = 0
    do k = poten%nz-1, 2, -1
        icount=0
        Pnot(k) = 0.d0 !Not counts from zero
        Icur(k) = Icur(k+1) / dS !upwards current is accumulative
        Ar(k) = 1.d-20
        point%kT = kBoltz * heat%tempinit(k) 
        do j=2,poten%ny-1
            do i=2,poten%nx-1
                ptype = classify_point(i,j,k)
                if (ptype == 0) then !surface point
                    point%Nstart = [i,j,k]
                    call emit_atpoint(poten,point)
                    if (point%ierr /= 0) then !something wrong in J coming out
                        if (debug) &
                            print *, 'WARNING: NaN or Inf found in emission calc'
                        Icur(k) = Icur(k) + 1.d-200
                        Nerrors = Nerrors + 1
                    else
                        Pnot(k) = Pnot(k) + point%heatNot !add Not heat 
                        Icur(k) = Icur(k) + point%Jem  ! add Current density to slice
                        Ar(k) = Ar(k) + dS
                        Ncorrect = Ncorrect + 1
                    endif
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
        if ((.not. intip) .and. (Icur(k) >0.d0)) then
            heat%tipbounds(2) = k
            intip = .true.
        endif
        if (intip) then
            rho = resistivity(heat%tempinit(k)) * fse
            Icur(k) = Icur(k) * dS !multiply by elementary area to get current
            Pnot(k) = Pnot(k) * dS !multiply by elementary area to get Watts
            Pjoule(k) =  rho * heat%dx * (Icur(k)**2) / max(Ar(k), 1.d-100)
            if (compmode >= 3) Pnot(k) = 0.d0 
            heat%hpower(k) =  (Pjoule(k) + Pnot(k)) / (Ar(k) * heat%dx)
            !Total heat density in SI units
        endif
            
        if (intip .and. Icur(k) <= 0.d0) then
            heat%tipbounds(1) = k + 1
            intip = .false.
        endif
    enddo
    heat%J_avg = Icur / Ar

    if (timings) print *, point%timings
    
    if (printheat) then
        print '(A15,F13.5,A10)', 'TipHeight =', &
            (heat%tipbounds(2) - heat%tipbounds(1)) * poten%grid_spacing(3), 'nm' 
        print '(A15,ES13.5,A10/A15,ES13.5,A10)', 'Icurbase =', &
                Icur(heat%tipbounds(1)),  'A', &
                'Pnot =', sum(Pnot(heat%tipbounds(1) : heat%tipbounds(2))), 'W'
        print '(A15,ES13.5,A10/A15,ES13.5,A10)', 'Pjoule =', &
            sum(Pjoule(heat%tipbounds(1) : heat%tipbounds(2))), 'W', &
            'Ptot@top =', heat%hpower(heat%tipbounds(2)), 'W/nm^3'
        print '(A15,F13.5,A10)', 'Fmax =', point%Fmax, 'V/nm'
        print '(A15,I12,A1,I10)', 'Nerr / Ncorr =', Nerrors , '/', Ncorrect
    endif
    
    if (savedata .and. mod(Nsteps,Nsavedata) == 0) then
   
        if(Ncalls == 0) then
            open(fid,file='out/heatdata.dat', action='write', status='replace')
            write (fid,*) 'Ncalls, Pjoule, Pnot, tempfinit, zi'
        endif
        
        write (fid,*)
        write(fid,*) Ncalls
        write (fid,*) Pjoule(heat%tipbounds(1):heat%tipbounds(2))
        write (fid,*) Pnot(heat%tipbounds(1):heat%tipbounds(2))
        write (fid,*) heat%tempinit(heat%tipbounds(1):heat%tipbounds(2))
        write (fid,*) &
            [(poten%grid_spacing(3)*i, i=1,(1+heat%tipbounds(2)-heat%tipbounds(1)))]
    endif

    Ncalls = Ncalls + 1
    
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
    use GeTElEC, only: EmissionData, cur_dens, print_data, plot_barrier
    
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
    logical                 :: badcontition, firsttime = .true.
    real(dp)                :: t1, t2

    badcontition = .false.
    if (firsttime) then
        allocate(that%xr(Nr), that%Vr(Nr))
        firsttime = .false.
    endif
    
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
    if (Fnorm >= point%Fmax * Flimratio .and. (compmode == 1 .or. compmode==3)) then
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
            
            if (Vmax > 1.1d0 * point%W) exit
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
        that%mode = -21 !in general case use polyfit trial
        if (badcontition) that%mode = -1 !use barrier model and force KX 
    else
        that%F = norm2(point%F)
        that%R = 1.d4
        that%gamma = 1.1d0
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
    
    if (that%ierr /= 0) then
        if (debug) then
            print *, 'Error in getelec cur_dens 1st call. Printing data:'
            print *, 'point:', point%Nstart
            call print_data(that, .true.)
        endif
        that%Jem = 0.d0
        that%heat = 0.d0
        point%ierr = 10+that%ierr
    endif
    
    if (that%Jem > point%Jmax * Jlimratio .and. (compmode==1 .or. compmode==3)) then
        that%full = .true.
        call cur_dens(that)
    
        if (that%ierr /= 0) then
            if (debug) then
                print *, 'Error in getelec cur_dens 2nd call. Printing data:'
                print *, 'point:', point%Nstart
                call print_data(that, .true.)
            endif
            that%Jem = 0.d0
            that%heat = 0.d0
            point%ierr = 20+that%ierr
        endif
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

    real(dp), parameter     :: Cv = 3.45d-6, L = 2.44d-8 
    ! Cv in(W*fs/(nm^3 K)), L=Wiedemann-Frantz constant in W*Ohm/K^2
    real(dp)                :: K, K1, rho 
    !K: Thermal conductivity, W/(nm K)
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

    do j = 1, solvesteps
        Told = Tnew
        do i = 2, Nx - 1
            rho = resistivity(Told(i)) * fse
            K = L * Told(i) / rho
            K1 = K / Cv
            a(i) = -K1 * heat%dt / (heat%dx**2)
            b(i) = 1.d0 + 2.d0 * K1 * heat%dt / heat%dx**2
            c(i) = -K1 * heat%dt / (heat%dx**2)
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
    heat%tempfinal(heat%tipbounds(2)+1:) = -1.d0 
  
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


subroutine poten_create(this, phi_in, gridspace)
    use bspline, only: db3ink

    type(Potential), intent(inout)  :: this
    real(dp), intent(in)            :: phi_in(:,:,:), gridspace(3)
    real(dp)                        :: x(size(phi_in, 1)), y(size(phi_in, 2)), &
                                        z(size(phi_in, 3))
    integer                         :: i, iflag
    real(dp)                        :: t2, t1 !timing
    
    if (debug) call cpu_time(t1)
    
    this%grid_spacing = gridspace
    if (this%nx/=size(x) .or. this%ny/=size(y) .or. this%nz/=size(z)) then
        this%nx = size(phi_in,1)
        this%ny = size(phi_in,2)
        this%nz = size(phi_in,3)
        if (allocated(this%phi)) deallocate(this%phi, this%bcoef, this%tx, &
                                            this%ty, this%tz)
        allocate(this%phi(this%nx, this%ny, this%nz), &
            this%bcoef(this%nx, this%ny, this%nz), &
            this%tx(this%nx + kx), this%ty(this%ny + ky), this%tz(this%nz + kz))
    endif
    x = [(this%grid_spacing(1)*(i-1), i=1,this%nx)]
    y = [(this%grid_spacing(2)*(i-1), i=1,this%ny)]
    z = [(this%grid_spacing(3)*(i-1), i=1,this%nz)]
    this%phi = phi_in
    call db3ink(x, this%nx, y, this%ny, z, this%nz, this%phi, kx, ky, kz, iknot, &
                this%tx, this%ty, this%tz, this%bcoef, iflag)
    if (iflag /= 0) then
        print *, 'Something went wrong with interpolation set. iflag=', iflag
        stop 
    endif
    this%set = .true.
    
    if (debug) then
        call cpu_time(t2)
        this%timing = this%timing + t2 - t1
    endif
end subroutine poten_create


function resistivity(T) result(rho)
    real(dp), intent(in)    :: T !temperature
    real(dp)                :: rho
    integer :: Ti
    include "rho_table.h"

    Ti = nint(T)/10

    if (Ti < 5 .or. Ti > 119) then
        if (debug) print *, &
            "WARNING: Temperature out of range in resisitivity calculation. T =", T

        Ti = min(Ti,119) ! cap up to 119
        Ti = max(Ti,5)   ! bottom down to 5
    end if
    ! Resistivity is an exponential function, so use log to make it linear, 
    !do normal interpolation, and then make it exponential again
    rho = exp( (log(rho_table(Ti+1)) - log(rho_table(Ti))) / (log(Ti*10d0+10)&
                - log(Ti*10d0)) * (log(T) - log(Ti*10d0)) + log(rho_table(Ti)) )

    rho = rho*10.d0 ! Convert to units Ohm*nm
end function resistivity


subroutine scale_vels(atoms, x1, heat, deltat_fs, delta, box)

    ! rewritten by A. Kyritsakis. First version not including Electromigration
    implicit none

    type(HeatData), intent(in)  :: heat
    
    integer, intent(in)     :: atoms(:,:,:) !Atoms on the grid 
    real(dp), intent(inout) :: x1(:) !< Atom velocities (parcas units)
    real(dp), intent(in)    :: deltat_fs, delta, box(3) 
    !Timestep (fs), Timestep inparcas units, Box size (Angstrom)


    real(dp)                :: Ek, Ek2, lambdas, sumlambdas
    !Ek: Kinetic energy a each k-slice
    !Ek0: Kinetic energy corresponding to temperature at each slice
    !lambdas: scaling factor for each k-slice
    !sumlambdas: sum of lambdas to print the average over tip
    integer             :: ilambdas ! counter of lambdas calcated for ti
    integer             :: i, j, k, m, Nslice
    !m: atom index on the xq velocities vector of parcas
    !Nslice : No of atoms in each slice
    
    sumlambdas = 0.d0
    ilambdas = 0

    do k = heat%tipbounds(1), heat%tipbounds(2)
        Ek = 0.d0 !put initial kinetic to zero for each slice
        Nslice = 0
        do j = 1, size(atoms,2)
            do i = 1, size(atoms,1)
                m = atoms(i, j, k)
                if (m > 0) then !if grid point has atoms
                    m = (m - 1)*3
                    Ek = Ek + ( x1(m + 1)**2 * box(1)**2 + &
                        x1(m + 2)**2 * box(2)**2 + &
                        x1(m + 3)**2 * box(3)**2) / (2.d0 * delta**2)
                    ! add kinetic energy of atom to total slice kinetic energy
                    Nslice = Nslice + 1
                endif
            enddo
        enddo
        Ek2 = 1.5d0 * Nslice * kBoltz * heat%tempfinal(k) 
        !Energy corresponding to slice temperature (Desired energy)
        lambdas = sqrt(1.d0 + deltat_fs * (Ek2 / Ek - 1.d0) / tau )
        !Scaling factor according to Brendsen theory
            
        if (Ek > 0) then ! if there were any atoms on the slice
            do j = 1, size(atoms, 2) !loop over all atoms of the slice
                do i = 1, size(atoms, 1)
                    m = atoms(i, j, k)
                    if (m > 0) then
                        m = (m - 1)*3
                        x1(m + 1 : m + 3) = x1(m + 1 : m + 3) * lambdas
                    endif
                    !scale velocities
                enddo
            enddo
            sumlambdas = sumlambdas + lambdas !Add sumlambdas to keep track
            ilambdas = ilambdas + 1
        end if
    enddo

    if(debug) print '(A15,ES13.5)', '<lambas> =', sumlambdas / ilambdas

end subroutine scale_vels

subroutine read_params()

    integer, parameter  :: fid = 549687
    character(len=12)   :: str
    
    open(fid, file = 'in/heatparams.in', action = 'read')
    
    read(fid,'(2A)')
    read(fid,*) str, workfunc
    if (str(1:8)/='workfunc') stop 'ERROR: not correct order in heatparams'
    read(fid,*) str, tau
    if (str(1:3)/='tau') stop 'ERROR: not correct order in heatparams'
    read(fid,*) str, compmode
    if (str(1:8)/='compmode') stop 'ERROR: not correct order in heatparams'
    close(fid)
    
    print *, 'heatparams read successfully', workfunc, tau, compmode

end subroutine read_params

end module heating
