module GeTElEC
!************************************************************************************
! General Tool for Electron Emission Calculations
!Author: Andreas Kyritsakis, University of Helsinki, 2016
!This module aims to calculate field emission current density and Nottingham effect 
!heating. The input has to be given in the EmissionData type structure. The output
!is also collected from there. All the theory behind this module can be found in 
!A. Kyritsakis, F. Djurabekova, https://arxiv.org/pdf/1609.02364.pdf
!Background theory can be found in
!A. Kyritsakis, J. Xanthakis, J. Appl. Phys. 119, 045303 (2016); 
!http://dx.doi.org/10.1063/1.4940721
! 
!************************************************************************************

use std_mat, only: diff2, local_min, linspace

implicit none
                                                                    
private
public  ::  cur_dens, C_wrapper, print_data, EmissionData, plot_barrier, debug, &
            dp, kBoltz, gamow_KX, gamow_num, gamow_general, cur_dens_SC, theta_SC

!************************************************************************************
!Global parameters not to be defined by the user

integer, parameter      :: dp = 8, sp = 4, fiderr = 987465, iknot = 0, &
                            idx = 0, knotx = 4 , fidparams = 812327
!Ny: length of special functions array
!knotx: No of bspline knots.
!idx, iknot: spline module parameters to bee kept 0
!dp, sp: double and single precision parameters
!fiderr: file id number for outputing error file
!above==0, the Transmission coefficient above Um is considered 1, above==1 Kemble 
real(dp), parameter     :: pi = acos(-1.d0), b = 6.83089d0, zs = 1.6183d-4, &
                           gg = 10.246d0, Q = 0.35999d0, kBoltz = 8.6173324d-5
!b, gg: exponential constants for Gamow
!kBoltz: Boltzmann constant
!zs : sommerfeld constant
! Q image potential constant
! all universal constants are in units nm, eV, Α


integer, parameter     :: SHARP = 1, BLUNT = 0, FIELD = 1, INTER = 0, THERMAL = -1

character(len=20), parameter   :: errorfile = trim('GetelecErr.txt'), &
                                  paramfile = trim('in/GetelecPar.in')
! names for the error output file and the parameters input file

!************************************************************************************
!Global parameters that are defined by the user from the params.txt file 

real(dp), save          :: xlim = 0.1d0, gammalim = 1.d3,  varlim = 1.d-3, &
                           epsfit = 1.d-4, nlimfield = 0.6d0, &
                           nlimthermal = 2.5d0, nmaxlim = 2.5d0                           
integer, save           :: Nmaxpoly = 10, debug = 0, above = 1
logical, save           :: spectra= .false., firstcall = .false.
character(len=50), save :: outfolder = "."

!xlim: limit that determines distinguishing between sharp regime (Numerical integral)
!and blunt regime where KX approximation is used
!nlimfield, nlimthermal are the limits for n to distinguish regimes
!nmaxlim: the maximum acceptable Jensen's n, above which the MG version of the 
!temperature corrected FN formula is used

!gammalim : maximum acceptable gamma. above it KX is forced

!epsfit : accuracy required in fittings
!varlim : limit of variance for the fitting approximation (has meaning for mode==-2)
!if the fitting variance is more than varlim we switch to spline mode

!Nmxpoly: max degree of the fitted polynomial
!spectra: set to true if you want to output spectroscopy data 
!if debug == 0, everything is silent. debug==1, error file is printed, 
! debug == 2 also some warning printed, parts are timed and calls are counted,
!if debug == 3 all warnings are printed and barrier is plotted 

!************************************************************************************
logical, save          :: readparams = .false.

type, public    :: EmissionData
!this type holds all the crucial data for the calculation of emission

    real(dp)    :: F=5.d0, R=100.d0, gamma=10.d0
        !Electrostatics: Field, Radius, enhancement factor
    real(dp)    :: W=4.5d0, kT=2.5d-2
        !Materical Characteristics: Work funciton,Temperature
    real(dp)    :: Jem=1.d-200, heat=1.d-200, Gam=1.d10 
        !Calculated results: Current density and Nott. heat
    real(dp)    :: theta = 1.d0
        !Space Charge field reduction factor
    real(dp)    :: xm=-1.d20, Um=-1.d20, maxbeta=0.d0, minbeta=0.d0, xmax = 2.d0, &
                    barlength = 1.d0
        !Barrier characteristics: xm=x point where barrier is max,
        ! Um=maximum barrier value, maxbeta=dG/dE @Fermi, minbeta=dG/dE @Um
        !xmax: estimation for maximum extent of the barrier
        !barlength: accurate barrier length between roots
    integer   :: regime = 1, sharpness = 0
        !1 for field, 0 for intermediate, -1 for thermal regimes
        !1 for sharp tip (numerical integration) and 0 for blunt (KX approx)
    
    real(dp), allocatable   :: xr(:), Vr(:)     !xr(nx), Vr(nx)
        !electrostatic potential externally defined as Vr(xr)
        
        
    real(dp), allocatable   :: tx(:), bcoef(:)  !tx(nx + kx), bcoef(nx)
        !bspline related parameters

    real(dp), allocatable   :: Apoly(:)
    integer                 :: Ndeg !degree of polynomial fit
        ! polynomial fitting working and result array created by dpolfit
    
    integer                 :: mode = 0
        !Mode of barrier calculation:  
        !0 : Barrier model (F,R,gamma)
        
        !1: Interpolation of Vr(xr). 
        !2: Interpolation and db1ink has been set up.
        
        
        !-1 : Barrier model, but force 'Blunt KX approximation'. "Bad" barrier.
    
        !-10: Fitting external data to (F,R,gamma). If unsatisfactory turn to FN
        !-11 : Fit to (F,R,gamma). If unsatisfactory, use interpolation
        !-12 : Fit to (F,R,gamma). Fitting is already done successfully
        
        !-20 : Fit to polynomial. If unsatisfactory turn to FN.
        !-21 : Fit to polynomial. If fitting not satisfactory, use interpolation.
        !-22 : Fit to polynomial and fitting has already been done, Apoly is ready.
        
    integer                 :: approx = 1
    !2: full calculation, force full integration. No GTF approximation for any regime
    !1: full calculation
    !0: GTF approximation
    !-1: FN approximation (Murphy - Good version)
    !-2: RLD approximation
        
    integer                 :: ierr = 0
        ! Integer to store error information.
        ! -1: invalid kT or W. Return zeroes without calculation 
        ! 0: Everything is fine and well - defined.
        ! 1: Mode should give xr, Vr and they are not allocated properly
        ! 2: Wrong input values for xr, Vr are given: non-monotonous
        ! 3: KX approximation giving negative U''(x)
        ! 4: Unknown error: some NaN appeared or Jem<0
        ! 5: NaN values over the energies more than 20%
        ! 10+: Some error in dfzero 1st appeared. + gives the ierr of dfzero
        ! 20+: Same as previous but for dfzero 2nd
        ! 30+: Save as previous but for dqage
        ! -2: Model fitting failed. Rough FN approximation  
        
    real(dp)                :: timings(5) = 0.d0
    !timing variables for cpu cost profiling
    !1: fitting time. 2: interpolation set time. 
    !3: J_num_integ time with barrier model. 4: J_num_integ with interpolation
    !5: GTF time.
    integer                 :: counts(5)  = 0
end type EmissionData

!************************************************************************************
contains


subroutine cur_dens(this)
!Calculates current density, main module function

    type(EmissionData), intent(inout)      :: this !main data handling structure
    
    real(dp)        :: Jf, Jt, n, s, E0, dE, heatf, heatt, F2, Fend, var
    real(dp)        :: nlimf, nlimt   !limits for n to distinguish regimes
    real(dp)        :: t1,t2,t3 !timing variables
    real(dp), allocatable :: xtemp(:), Vtemp(:)
    integer         :: i, GTFerr
    character(len = 50) :: msgerr
    
    if (debug > 1) then
        call cpu_time(t1)
        print *, 'entering cur_dens'
    endif
    
    if (.not. readparams) call read_params()
    
    this%ierr = 0
    GTFerr = 0
    
    if (this%W <= 0.d0 .or. this%kT < 0 .or. isnan(this%W) .or. isnan(this%kT)) then 
        !check input validity
        this%ierr = -1
        this%Jem = 0.d0
        this%heat = 0.d0
        call error_msg(this,'Wrong input data.')
        return
    endif
    
    !preparing calculation according to calculation mode
    if (this%mode > 0 .or. this%mode < -1) then
        if (size(this%xr) < 2 .or. size(this%Vr) /= size(this%xr) .or. & 
                    .not. (allocated(this%xr) .and. allocated(this%Vr))) then
            call error_msg(this,'xr,Vr not properly allocated')
            this%ierr = 1
            return
        else !everything ok with sizes of xr, Vr
            do i = 2,size(this%xr) !iterate over xr,Vr to check if monotonous
                if (this%xr(i) < this%xr(i-1) .or. this%Vr(i) < this%Vr(i-1)) then
                    if (this%Vr(i) > this%W) then !if length enough to cover barrier
                        allocate(xtemp(i-1), Vtemp(i-1))
                        xtemp = this%xr(1:i-1)
                        Vtemp = this%Vr(1:i-1)
                        deallocate(this%xr, this%Vr)!get rid of the last bad values
                        allocate(this%xr(i-1), this%Vr(i-1))
                        this%xr = xtemp
                        this%Vr = Vtemp
                        this%ierr = 0 !return to normal execution
                        exit
                    else
                        this%ierr = 2
                        write(msgerr, &
                            '("Non-monotonous xr,Vr at i=",i4,". Not recovered")') i
                        call error_msg(this, trim(msgerr))
                        return
                    endif
                endif
            enddo
        endif
    endif
    
      
    if (this%mode < -1) then !fit external data
        if (this%mode == -10 .or. this%mode == -11) then
            var = fitpot(this) ! fit to (F,R,gamma) model
        elseif (this%mode == -20 .or. this%mode == -21) then
            var = fitpoly(this) ! fit to polynomial
        endif
        
        if ( var < varlim .and. (.not. isnan(var))) then
            this%mode = (this%mode / 10) * 10 - 2 ! set to -22 or -12
        elseif (mod(this%mode, 10) == -1 .and. this%Vr(size(this%Vr)) > this%W) then
            this%mode = 1 !switch to interpolation
            if (debug > 1) &
                print *, 'var =', var,'varlim =',varlim, 'Switched to interpolation.'
        else ! mode = -10/20 or Vr(xr) not covering barrier, turn to error catching
            this%mode = -1
            this%ierr = -2
            this%R = 1.d4
            this%gamma = 1.d0
            call error_msg(this,'Warning: Fitting failed. Rough FN approximation') 
        endif
    endif
    
    
    if (this%mode == 1 .or. this%mode == -22) then        
        this%F = ( this%Vr(2) - this%Vr(1) ) / ( this%xr(2) - this%xr(1) )
        F2 = ( this%Vr(3) - this%Vr(2) ) / ( this%xr(3) - this%xr(2) )
        Fend = ( this%Vr(size(this%Vr)) - this%Vr(size(this%Vr)-3) ) / &
                ( this%xr(size(this%xr)) - this%xr(size(this%xr)-3) )
        this%R= abs(this%F/((F2-this%F)/(this%xr(3)-this%xr(1))))
        this%gamma = this%F / Fend
        
        if (debug > 1) then
            print *, "Parameters were guessed:"
            call print_data(this)
        endif
    endif
    
    if (debug > 1) then !calculate preparation time
        call cpu_time(t2)
        this%timings(1) = this%timings(1) + t2 -t1
        this%counts(1) = this%counts(1) + 1
    endif

    if (this%gamma < 0.d0) then !check if gamma is unacceptable
        this%gamma = 1.1d0
        call error_msg(this, 'Negative gamma found. Error recovered with gamma=1.1')
    elseif (this%gamma > gammalim) then 
        this%gamma = gammalim
        call error_msg(this, 'Too high gamma found. Recovered with gamma=gammalim')
    endif
    
    if (this%R < 0.1d0) then
        if (debug > 1) print *, 'R < 0.1d0 was inputted. Recovered with R = 0.1'
        this%R = 0.1d0
    endif
    
    if (this%approx == 1) then !if full calculation expand borders of regimes
        nlimf = nlimfield
        nlimt = nlimthermal
    elseif (this%approx == 2) then !force full integration always
        nlimf = -100.d0
        nlimt = 1.d100
    else !if only GTF equation standard limit = 1
        nlimf = 1.d0
        nlimt = 1.d0
    endif
    
    this%Um = -1.d20  !set maximum values to unknown
    this%xm = -1.d20
    
    call gamow_general(this,.true.) !calculate barrier parameters
    
    if (debug > 2.) then
        call plot_barrier(this)
    endif
    
    if (this%ierr > 0) then 
        call error_msg(this,'Error after first gamow_general.')
        return
    endif
    
    if (this%approx >= 0) then !full calculation or GTF approximation
        if (this%kT * this%maxbeta < nlimf .and. &
                (this%Um > 0.5 .or. (this%approx == 0))) then!field regime
            n = 1.d0/(this%kT * this%maxbeta) !standard Jensen parameters
            s = this%Gam
            this%regime = FIELD
            
        else if ((this%kT * this%minbeta > nlimt .or. this%maxbeta > 1.d10 ).and. &
                    (this%Um > 0.5 .or. (this%approx == 0))) then !thermal regime
            n = 1.d0/(this%kT * this%minbeta) !standard Jensen parameters
            s = this%minbeta * this%Um
            this%regime = THERMAL
        else  !intermediate regime
            this%regime = INTER
            if (this%approx >= 1) then
				call J_num_integ(this) !Numerical integration over energies
            else
                GTFerr = GTFinter(this)
            endif
            if (GTFerr/=0) call J_num_integ(this)  !if GTF return with error, do full
        endif
        
        if (debug > 1) print '(A10, ES12.4/A10, ES12.4/A10, ES12.4/A10, ES12.4)', &
                'n =', n, 's =', s, 'Sig(n) =', Sigma(n), 'Sig(1/n) =', Sigma(1/n) 
        
        if (this%regime /= INTER) then
            if (n > nmaxlim) then !for n>5 use MG expressions. Approximations for 
                               ! Sigma misbehave  
                this%Jem = zs*pi* this%kT * exp(-this%Gam) & !Murphy-Good version FN
                    / (this%maxbeta * sin(pi * this%maxbeta * this%kT))
                this%heat = -zs*(pi**2) * (this%kT**2) * & !Nottingham at F regime
                    exp(-this%Gam)*cos(pi* this%maxbeta * this%kT)  &
                    / (this%maxbeta * (sin(pi * this%maxbeta * this%kT )**2))
            else !
                this%heat = zs * (this%kT**3) * ( ( (n*s + 1.d0) * Sigma(n) &
                    - DSigma(n)) * exp(-n*s) + (n**3) * DSigma(1.d0/n) * exp(-s))
                this%Jem = zs * (this%kT**2) * ((n**2) * Sigma(1.d0/n) * exp(-s) &
                    + exp(-n*s) * Sigma(n))
            endif
        endif
    elseif (this%approx == -1) then !Murphy-Good version FN
        this%Jem = zs * (this%maxbeta**(-2)) * exp(-this%Gam) * &
                (1.+ (pi * this%maxbeta * this%kT)**2 / 6.)
        this%heat = -zs * (this%maxbeta**(-3)) * exp(-this%Gam) * &
                (1. - (pi * this%maxbeta * this%kT)**2 / 6.)
    elseif (this%approx == -2) then ! RLD equation
        this%Jem = zs * (this%kT**2) * exp(-this%Um/this%kT) * &
                (1. + (pi * this%minbeta * this%kT)**(-2) / 6.) 
        this%heat = (this%kT + this%Um) * this%Jem
    else
        print *, 'Wrong approx value. this%approx = ', this%approx
        stop
    endif
    
    if (debug > 1) call print_data(this)
        
    
    this%heat = - this%heat !convention: consider deposited heat (+ for heating) 
    
    if ((isnan(this%Jem) .or. isnan(this%heat) .or. this%Jem < -1.d-201)  &
        .and. this%ierr == 0) this%ierr = 4
    
    if (this%ierr > 0) then 
        call error_msg(this,'Error in the end of execution.')
    endif
    
    contains

    pure function Sigma(x) result(Sig)!Jensen's sigma
        real(dp), intent(in)    :: x
        real(dp)                :: Sig
        
        if (x > 3.d0) then 
            Sig = 0.d0
        else
            Sig = (1.d0+x**2)/(1.d0-x**2)-.3551d0*x**2-0.1059d0*x**4-0.039*x**6
        endif
    end function Sigma
    
    pure function DSigma(x) result(ds) !Jensen's sigma'
        real(dp), intent(in):: x
        real(dp):: ds
        
        if (x > 3.d0) then
            ds = 0.d0
        else
            ds = (-1.d0 + 3.6449 * x**2 + 1.3925d0 * x**4 + 0.0853d0 * x**6 + &
                0.0723d0 * x**8 - 0.195d0 * x**10)/((-1.d0 + x**2)**2)
        endif
    end function DSigma

end subroutine cur_dens 


subroutine cur_dens_SC(this, Voltage)
!calculates the current density taking into account the space charge
! the planar diode model is used with voltage as external parameter
    type(EmissionData), intent(inout)       :: this !main data handling structure
    real(dp), intent(in)                    :: Voltage
    !voltage of the diode and previous guess for the theta factor
    
    integer     :: i
    real(dp)    :: theta_new, error
    
    
    do i = 1, 50
        if (debug > 0) print *, "SC cycle = ", i, "Voltage = ", Voltage
        
        !multiply multiply field by theta and calculate current density
        this%F = this%F * this%theta
        if (allocated(this%Vr)) this%Vr = this%Vr * this%theta
        call cur_dens(this)
        if (debug > 0) call print_data(this)
        
        theta_new = theta_SC(this%Jem, Voltage, this%F)
        error =  (theta_new - this%theta)
        if (debug > 0) print *, 'theta_new = ', theta_new, 'SC error = ', error
        
        !restore electrostatics
        this%F = this%F / this%theta
        if (allocated(this%Vr)) this%Vr = this%Vr / this%theta
        
        !calculate new reduction factor
 
        if (abs(error) < 1.e-5) return
        this%theta = 0.5 * (theta_new + this%theta) 
    enddo

end subroutine cur_dens_SC


subroutine gamow_general(this,full)
!Calculates Gamow exponent in the general case: Choose appropriate approximation
    use bspline, only: db1ink
    type(EmissionData), intent(inout)       :: this !data structure
    type(EmissionData)                      :: new
    logical, intent(in)                     :: full
    !determine if maxbeta, minbeta will be calculated. Not needed when numerical
    !integration over energies is done in the intermediate regime

    real(dp)                                :: dw!delta w.f. for maxbeta differential
    real(dp)                                :: x, xmaxallowed
    integer                                 :: iflag , info
    real(dp)                                :: t1, t2 !timing
    
    xmaxallowed = 1.d2
    x = this%W / (this%F * this%R)!length of barrier indicator
    dw = 1.d-2
    
    if (this%mode == 0 .or. this%mode == -12 .or. this%mode == -1) then
        if (x > 0.4d0 .and. this%gamma > 1.05d0) then !the second root is far away
            this%xmax = min(this%W * this%gamma / this%F, xmaxallowed)
            !maximum length xmaxallowed
        else  !the second root is close to the standard W/F
            this%xmax = 2.5d0 * this%W / this%F
        endif
    else
        this%xmax = this%xr(size(this%xr))
    endif
    
    
    if (x > xlim .and. this%mode /= -1) then !not x<<1, use numerical integration
        if (this%mode == 1) then    !setup spline interpolation
            if (debug > 1) call cpu_time(t1)
            if (.not. allocated(this%tx) .or. .not. allocated(this%bcoef) .or. &
                    size(this%tx) /= size(this%Vr) + knotx) then
                if (allocated(this%tx)) deallocate(this%tx)
                if (allocated(this%bcoef)) deallocate(this%bcoef)
                allocate(this%tx(size(this%Vr) + knotx), this%bcoef(size(this%Vr)))
            endif
            
            call db1ink(this%xr, size(this%xr), this%Vr, knotx, iknot, &
                        this%tx, this%bcoef, iflag)
            if (iflag/=0) print *, 'error in spline setup, iflag = ', iflag
            
            if (debug > 1) print *, 'spline interpolation is set'
            this%mode = 2 !change mode so no need to call db1ink again
            
            if (debug > 1) then
                call cpu_time(t2)
                this%timings(2) = this%timings(2) + t2 - t1
                this%counts(2) = this%counts(2) + 1
            endif 
        endif
        
        call gamow_num(this, .true.)
        this%sharpness = SHARP
        if (this%Um < 1.d-2 ) then!.or. (.not. full)) then !barrier lost
            this%maxbeta = this%minbeta
        elseif (this%Gam == 1.d20) then !too long barrier
            this%maxbeta = 1.d20
        elseif (full) then
            new = this !copy structure for calculation for E-dw
            new%W = this%W - dw
            call gamow_num(new,.false.)
            this%maxbeta = abs(new%Gam - this%Gam) / dw !derivative
        endif
    else !large radius, KX approximation usable
        
        info = gamow_KX(this,full)
        this%sharpness = BLUNT
        if (info == -2) this%ierr = 3
    endif
    
    if ((isnan(this%Gam) .or. isnan(new%Gam)) .and. debug > 1 ) then
        print *, 'Gam takes NaN in gamow_general. this ='
        call print_data(this)
        print *, 'new ='
        call print_data(new)
    endif
    
end subroutine gamow_general


subroutine gamow_num(this, full)
    !integrates numerically the barrier and obtains G
    ! if not full, Um, xm are already found and need not be calculated again
    use bspline, only: db1val
    use std_mat, only: binsearch

    type(EmissionData), intent(inout)   :: this
    logical, intent(in)                 :: full
        
    !integration and root finding tolerance parameters
    integer, parameter                  :: maxint = 256 !maximum integration intervs
    real(dp), parameter                 :: AtolRoot=1.d-6, RtolRoot=1.d-6, &
                                           AtolInt=1.d-5, RtolInt=1.d-5

    ! integration and root finding working variables, see slatec reference
    real(dp), dimension(maxint)         :: ALIST, BLIST, RLIST, ELIST
    real(dp)                            :: ABSERR
    integer                             :: IFLAG, NEVAL, IER, IORD(maxint), LAST
        
    real(dp)                            :: G, x = 1.d-5 , x1(2), x2(2), dx
    integer                             :: i, binout(2), indxm, ixrm !
                                        ! indxm: index of xr closest to xm
    
    ixrm = size(this%xr)
    
    if (full) then ! Um is not initialized
        this%Um = -local_min(x, this%xmax, 1.d-8, 1.d-8, neg_bar, this%xm)
        !Catch the case global minimum not found
        if (this%Um < 0.d0 .or. this%Um < bar(.5d0 * this%xm)) &
            this%Um = -local_min(x, this%xm-1.d-3, 1.d-8, 1.d-8, neg_bar, this%xm)

        if (this%mode <= 0) then  ! if model or polynomial
            dx = 1.d-2  !dx : the differentiation dx used for minbeta
        else !if interpolation choose big dx for diff2. Αvoid num instabillity
            binout = binsearch(this%xr,this%xm)  !
            if (binout(2) /= 0) then !'xr not sorted or xm out of bounds'
                this%ierr = 5
                return
            endif
            indxm = binout(1)
            if (indxm == 1) then  !avoid segfault
                dx  =  this%xr(2)-this%xr(1)
            else
                dx = this%xr(indxm) - this%xr(indxm -1)
            endif
        endif
        this%minbeta = 22.761d0 / sqrt(abs(diff2(bar,this%xm,dx)))
        if (isnan(this%minbeta) .and. debug > 1) print *, 'indxm=', indxm, &
                'dx=', dx, 'U(xm+dx)=', bar(this%xm+dx), 'U(xm-dx)=', bar(this%xm-dx)
        if (this%Um < 1.d-2) then !lost barrier
            this%Gam = this%minbeta * this%Um
            this%maxbeta = this%minbeta        
            return
        endif
        
                
    endif

    x1 = [0.01d0, this%xm] !interval for search of first root
    x2 = [this%xm, this%xmax] !interval for search of second root
    
    if (bar(this%xmax)  > 0.d0) then
        this%Gam = 1.d20
        return
    endif
    
    call dfzero(bar,x1(1),x1(2),x1(1),RtolRoot,AtolRoot,IFLAG) !first root
    if (IFLAG /= 1) then !if something went wrong
        this%ierr = 10 + IFLAG !set error flag
        if (debug > 1) &
            print *, 'Error occured in dfzero1. IFLAG =', IFLAG! print error messages
    endif
    
    call dfzero(bar,x2(1),x2(2),x2(1),RtolRoot,AtolRoot,IFLAG) !second root
    if (IFLAG /= 1) then !if something went wrong
        this%ierr = 20 + IFLAG !set error flag
        if (debug > 1)  &! print error messages
            print *, 'Error occured in dfzero2. IFLAG =', IFLAG
    endif
    
    if (this%ierr > 10 .and. debug  > 2) then
        print *, 'xm =', this%xm, '|| xmax =', this%xmax
        print *, 'x1=', x1, '|| x2=', x2
        print *, 'bar(x1)=', bar(x1(1)), bar(x1(2)), &
                    '|| bar(x2)=', bar(x2(1)), bar(x2(2))
    endif
    
    call dqage(sqrt_bar, maxval(x1), minval(x2), AtolInt,RtolInt, 2, maxint, &
                G, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
                !integrate
    if (IER /= 0) then
        this%ierr = 30 + IER
        if (debug > 2) then
            print *, 'Error occured in dqage. IER =', IER
            print *, 'limits = [', x1, '|| ', x2, ']'
            print *,  'Neval =', NEVAL, '| ABSERR =', ABSERR, '| LAST=', LAST 
        endif
    endif
    !caution in the interval limits. x1, x2 are not necessarily returned from dfzero
    !in increasing order. Use maxval, minval!!!
    
    this%Gam = gg * G
    
    this%barlength = x2(1) - x1(2)
    
    contains
    ! barrier functions needed to be passed as externals to dfzero and dqage
    
    function bar(x) result(V)!sphere barrier model
        real(dp), intent(in)    :: x
        real(dp)                :: V, Vinterp(2)
        integer                 :: iflag, inbvx
        
        select case (this%mode)
            case (2) !interpolation mode.. Use splines
                inbvx = 1
                call db1val(x, idx, this%tx, size(this%Vr), knotx, &
                        this%bcoef, Vinterp(1), iflag, inbvx)
                if (iflag /= 0) then !out of bounds.. do linear extrapolation 
                    Vinterp(1) = this%Vr(ixrm) + (x - this%xr(ixrm)) * &
                    (this%Vr(ixrm)-this%Vr(ixrm-1)) / (this%xr(ixrm)-this%xr(ixrm-1))
                endif
                V = this%W - Vinterp(1) - Q / (x + ( 0.5d0 * (x**2)) / this%R)
            case (-22) ! fitted polynomial mode. Evaluate fitted polynomial
                if (x<= this%xr(ixrm)) then
                    call dp1vlu(this%Ndeg, 0, x, Vinterp(1), Vinterp(2:), this%Apoly)
                else !out of bounds. Polynomial misbehaves. Use extrapolation
                    call dp1vlu(this%Ndeg, 1, this%xr(ixrm), & !polyval + der in xrm
                            Vinterp(1), Vinterp(2:), this%Apoly)
                    Vinterp(1) = Vinterp(1) + Vinterp(2) * (x - this%xr(ixrm))
                    !linear extrapolation using the derivative of the polynomial
                    Vinterp(1) = Vinterp(1)
                endif
                V = this%W - Vinterp(1) - Q / (x + ( 0.5d0 * (x**2)) / this%R)
            case (0,-12, -1)
            !use standard F,R,gamma model for the electrostatic potential
                V = this%W - (this%F * this%R * x*(this%gamma - 1.d0)  &
                + this%F * x**2) / (this%gamma * x + this%R * (this%gamma - 1.d0)) &
                - Q / (x + ( 0.5d0 * (x**2)) / this%R)
            case default
                print *, 'Error: invalid mode in bar(). Mode = ', this%mode  
                stop
        end select
    end function bar
    
    function neg_bar(x) result(nv)
        real(dp), intent(in)    ::x
        real(dp)                ::nv
        nv = - bar(x)
    end function neg_bar
    
    function sqrt_bar(x) result(rv)
        real(dp), intent(in)    ::x
        real(dp)                ::rv
        rv=sqrt(bar(x))
    end function sqrt_bar
    
end subroutine gamow_num


function gamow_KX(this, full) result(info) 
    !calculate Gamow parameters using Kyritsakis-Xanthakis approximation
    use std_mat, only: lininterp
    use ellfuns, only: Ny, vv, tt, ww, psi
    
    type(EmissionData), intent(inout)   :: this
    logical, intent(in)                 :: full !if F, only Gamow is calculated
    
    real(dp)                            :: yf, t, ps, v, omeg, temp, chi
    integer                             :: info
    !info : information integer showing different extreme cases
    !1: yf>1. Linear extrapolation used
    !-1: yf<0. Something is seriously wrong
    !0: everything is ok 0<yf<1
    !2: Um<0. G and maxbeta are calculated according to minbeta lin extrapolation
    !-2: U''(x) is negative. The approximation is invalid. Switch to numerical 

    yf = 1.44d0 * this%F / (this%W ** 2)
    chi = this%W / (this%F * this%R)
    if (yf>1.d0) then
        t = tt(Ny) + (yf - 1.d0) * (tt(Ny) - tt(Ny-1))*(Ny - 1.d0)
        ps = psi(Ny) + (yf - 1.d0) * (psi(Ny) - psi(Ny-1))*(Ny - 1.d0)
        v = vv(Ny)+ (yf - 1.d0) * (vv(Ny) - vv(Ny-1))*(Ny - 1.d0)
        omeg = ww(Ny)+ (yf - 1.d0) * (vv(Ny) - vv(Ny-1))*(Ny - 1.d0)
        info = 1
    elseif  (yf<0.d0) then
        t= tt(1)
        ps = psi(1)
        v = vv(1)
        omeg = ww(1)
        info = -1
    else
        t = lininterp(tt,0.d0,1.d0,yf)
        ps = lininterp(psi,0.d0,1.d0,yf)
        v = lininterp(vv,0.d0,1.d0,yf)
        omeg = lininterp(ww,0.d0,1.d0,yf)
        info = 0
    endif

    if (full) then
        this%Um = this%W - 2.d0*sqrt(this%F * Q) + 1.5d0 * Q / this%R
        this%xm = sqrt(Q / this%F) + Q / (this%F * this%R)
        temp =  (this%F**1.5d0) / sqrt(Q) - 4.d0 * this%F / this%R
        if (temp < 0.d0) then
            if(debug > 2) print *, 'Warning: minbeta goes to negative root'
            temp  =  1.d-10  !make sure that sqrt is positive, avoid NaNs
            info = -2
        endif
        this%minbeta = 16.093d0/sqrt(abs(temp))
        if (this%Um > 0.d0) then
            this%maxbeta = gg*(sqrt(this%W)/this%F) * &
                (t+ps*(this%W / (this%F * this%R)))
        else
            this%maxbeta = this%minbeta
        endif
    endif
    
    
    if (this%Um > 0.d0) then
        this%Gam = (2*gg/3)*((this%W **1.5d0) / this%F) &
            * (v+omeg*(this%W / (this%F * this%R)))
        this%barlength = (this%W / this%F)  * &
                    (sqrt(1-yf) + chi * (1.d0 - .625*yf)/sqrt(1-yf))
    else
        this%Gam = this%minbeta * this%Um
        info = 2
    endif

end function gamow_KX


function GTFinter(this) result(error)
!Calculates estimation of current according to GTF theory for intermediate regime.
!It is used when a fast calculation is needed and we don't care about full numerical
!calculation for all energies in the intermediate regime.

    type(EmissionData), intent(inout)   :: this
    
    ! calculated in single precision because slatec's polynomial root subroutine
    ! supports only sp
    real(sp)            ::  polybeta(3), work(8)
    complex(sp)         :: rt(2)
    integer             :: ierr, error
    real(dp)            :: zmax, Em, Gm, s, C_FN, Bq, B_FN, UmaxokT
    real(dp)            :: t1, t2 !timing variables
    
    error = 0    
    if (debug > 1) call cpu_time(t1)
    C_FN = this%maxbeta * this%Um
    Bq = this%minbeta * this%Um
    B_FN = this%Gam
    UmaxokT = this%Um / this%kT !Umax/kT
    
    polybeta = real([3.*(C_FN+Bq-2.*B_FN), 2.*(3.*B_FN-Bq-2.*C_FN), C_FN-UmaxokT],sp)
    !polynomial describing G(E)
    if (abs(polybeta(1)) <1.d-2 .and. debug > 1) then
        print *, 'error: GTF polynomial with negative root'
        call print_data(this)
    endif
    
    call RPQR79(2,polybeta,rt,ierr,work) !finding the roots of the polynomial
    if (ierr /= 0) then
        print *,'polynomial root finder in GTFinter encounter error. printing data:'
        call print_data(this)
    endif
    
    
    zmax=minval(real(rt,dp)) ! take the smaller root
    if (debug > 1) print *, 'Polynomial roots in GTFinter: rt = ', rt
    if (zmax<0) then
        zmax=maxval(real(rt,dp)) !choose the positive root if it is negative
    endif
    
    Em = this%Um * zmax !E where max integrant occurs according to Jensen
    Gm = -(C_FN+Bq-2.*B_FN) * zmax**3 - (3.*B_FN-Bq-2.*C_FN) * zmax**2 &
            -C_FN * zmax + B_FN ! G at the maximum
    if (Gm<0.d0 .or. Gm > B_FN) then
        if (debug > 1) print *, 'Error in GTFinter, Gm =', Gm
        error = 1
        return
    endif
    
    s = Em / this%kT + Gm
    if (Gm < 0.d0 .and. debug > 1) then
        print *, 'zmax = ', zmax, 'Gm = ', Gm
        call print_data(this)
    endif
    
    this%heat = zs*(this%kT**3)*exp(-s) * (s + 1.d0 + .5d0*s**2)
    this%Jem = zs * (this%kT**2) * exp(-s) * (s + 1.d0)
    if (debug > 1) then
        call cpu_time(t2)
        this%timings(5) = this%timings(5) + t2 -t1
        this%counts(5) = this%counts(5) + 1
    endif

end function GTFinter


subroutine J_num_integ(this)
!numerical integration over energies to obtain total current according
!to landauer formula
    type(EmissionData), intent(inout)   :: this
    type(EmissionData)                  :: new
    
    real(dp), parameter                 :: cutoff=1.d-4 
            !cutoff of the exponentially decreasing
    integer, parameter                  :: Nvals=512, fidout=1953
            !no of intervals for integration and fidout for spectroscopy
            
    real(dp)                            :: Gj, G(4*Nvals), Ej, dE, Ea, Eb
    real(dp)                            :: insum, Jsum, Jcur, Umax, Emax, Emin, E0
    real(dp)                            :: Gup(2*Nvals), Gdown(2*Nvals), fj, fmax
    real(dp)                            :: integrand, outsum, t2, t1
    integer                             :: j, i, k, NaNs

    if (debug > 1) call cpu_time(t1)
    if (spectra) then !if spectra then j(E) will be printed into a file
        open(fidout,file=trim(outfolder)//'/spectra.csv')
        write(fidout,*) '# Energy [eV],    Nott. spectral dens. [W /nm^2 /eV],', &
            '  Spectral cur. dens. [A/nm^2 /eV],', ' Fermi-Dirac integral log,', ' Gamow exponent'
    endif
    this%Um = -1.d20
    this%xm = -1.d20
    
    call gamow_general(this,.true.) !make first estimation and calculate Um, xm etc.
    
    if (this%Gam == 0.d0) then!if no barrier for fermi electrons
        this%Jem = 1.d20
        return
    endif

    if (this%kT * this%maxbeta < 1.05d0) then!field regime, max is close to Ef
        Emax = abs(log(cutoff)) / (1.d0 / this%kT - .7 * this%maxbeta)
        Emin = - abs(log(cutoff)) / this%maxbeta
        E0 = 0.d0
    else if (this%kT * this%minbeta > .95d0) then!thermal regime, max is close to Um
        Emax = this%Um + abs(log(cutoff)) * this%kT
        Emin = this%Um - abs(log(cutoff)) / (this%minbeta - 0.7 / this%kT)
        E0 = this%Um
    else! integmediate regime, max is in between ,find it
        Ea = 0.d0
        Eb = this%Um
        do j=1,100!bisection to find where n=1
            Ej=(Eb+Ea)*.5d0
            Umax = this%Um - Ej
            new = this
            new%W = this%W - Ej
            call gamow_general(new,.true.)
            if (this%kT * new%maxbeta < 0.98d0) then
                Eb = Ej
            else if (this%kT * new%maxbeta > 1.02d0) then
                Ea = Ej
            else
                exit
            endif
        enddo
        Emax = this%Um
        Emin = 0.d0
        E0 = Ej
    endif
    
    dE = (Emax-Emin) / (Nvals-1)
    Jsum = 0.d0
    insum = 0.d0
    outsum = 0.d0
    fmax = 0.d0
    NaNs = 0
    
    do j = 1, 2*Nvals !integrate over high energies
        Ej = E0 + (j-1) * dE
        Umax = this%Um - Ej
        if (Umax > 0.d0) then
            new = this
            new%Um = Umax
            new%W = this%W - Ej
            call gamow_general(new,.false.) !fastly calculate only G 
            Gj = new%Gam
        elseif (above /= 0) then
            Gj = this%minbeta * Umax
        else
            Gj = -200.d0
        endif
        if (isnan(Gj)) then
            NaNs = NaNs + 1
            if (debug > 1) print *, 'warning: NaN in Gj for Ej  = ', Ej
            if (debug > 2) call error_msg(this,'NaN in Gj')
            if (j==1) then
                Gup = this%Gam
            else
                Gup(j) = Gup(j-1)
            endif
        else
            Gup(j) = Gj
        endif
        fj=lFD(Ej,this%kT)/(1.d0+exp(Gup(j)))
        if (fj > fmax) fmax = fj
        if (abs(fj/fmax) < cutoff) exit
        Jsum = Jsum + fj
    enddo
    if (j>2*Nvals) j = 2*Nvals
    do i=1,2*Nvals!integrate over low energies
        Ej = E0 - i * dE
        Umax = this%Um - Ej
        if (Umax>0.d0) then
            new = this
            new%W = this%W - Ej
            new%Um = Umax
            call gamow_general(new,.false.)
            Gj = new%Gam
        elseif (above /= 0) then
            Gj = this%minbeta * Umax
        else
            Gj = -200.d0
        endif
        if (isnan(Gj)) then
            NaNs = NaNs + 1
            if (debug > 1) print *, 'warning: NaN in Gj for Ej  = ', Ej
            if (debug > 2) call error_msg(this,'NaN in Gj')
            if (i==1) then
                Gdown = this%Gam
            else
                Gdown(i) = Gup(i-1)
            endif
        else
            Gdown(i) = Gj
        endif
        fj=lFD(Ej,this%kT)/(1.d0+exp(Gdown(i)))
        if (fj>fmax) fmax=fj
        if (abs(fj/fmax) < cutoff) exit
        Jsum = Jsum + fj
    enddo
    
    if (NaNs > (i+j) / 5) then ! if too many NaNs appeared, the result is wrong
        this%ierr = 5
        this%heat = 0.d0
        this%Jem = 0.d0
        return
    endif
    
    this%Jem = Jsum * zs * this%kT * dE
    if (i>2*Nvals) i = 2*Nvals

    G(1:i+j) = [Gdown(i:1:-1),Gup(1:j)]!bring G vectors together and in correct order

    do k=1,i+j!integrate for nottingham effect
        if (above == 0 .and. Ej > this%Um) G(k) = -100.d0
        integrand = Ej * (insum + .5d0*dE / (1.d0 + exp(G(k)))) & 
                    / (1.d0 + exp(Ej/this%kT))
        !.5d0 ... term is what remains from insum for ending trapezoid rule
        insum = insum + dE / (1.d0 + exp(G(k)))
        outsum = outsum + integrand
        if (spectra) then !write spectroscopy data into the spectroscopy file
            fj=lFD(Ej,this%kT)/(1.d0+exp(G(k)))
            write(fidout,*) Ej, zs * integrand, zs * integrand/Ej, fj, G(k)
        endif
        Ej = Ej + dE
    enddo

    this%heat = zs * outsum * dE
    
    if (spectra) then
        close(fidout)
    endif
    
    if (debug > 1) then
        call cpu_time(t2)
        i = 3
        if (this%mode > 0) i = 4 !check if interpolation or barrier model
        this%timings(i) = this%timings(i) + t2 - t1
        this%counts(i) = this%counts(i) + 1
    endif

    contains
    
    pure function lFD(E,kT) result(L) !Fermi dirac log function 
        real(dp), intent(in)::E,kT
        real(dp) :: L
        if (E>10.d0*kT) then
            L=exp(-E/kT)
        else if (E < -10.d0 * kT)  then
            L = - E / kT
        else
            L=log(1.d0+exp(-E/kT))
        endif
    end function lFD
    
end subroutine J_num_integ


subroutine print_data(this, full, filenum)
    !print the state of the object nicely
    type(EmissionData), intent(inout)   :: this
    integer, intent(in), optional       :: filenum
    logical, intent(in), optional       :: full
    integer                             :: fid, i
    
    character(len=12)   :: reg_str(-1:2) = ['Thermal     ', &
                                            'Intermediate', &
                                            'Field       ', &
                                            'Unknown     ']
                                        
    character(len=10)    :: sharp_str(0:2) = [  'Blunt', &
                                                'Sharp', &
                                                'N/A  ']
                                            
    character(len=20)   :: approx_str(-2:2) = [ 'Richardson-Dushman  ', &
                                                'Murphy-Good         ', &
                                                'General T-F (Jensen)', &
                                                'Automatic selection ', &
                                                'Full integration    ' ]
                                                
    if (this%regime > 1 .or. this%regime < -1) this%regime = 2
    if (this%sharpness > 1 .or. this%sharpness < 0) this%sharpness = 2
    if (this%approx > 2 .or. this%approx < -2) stop 'Approximation set wrongly'
    
    if (.not. present(filenum)) then
        fid = 6
    else
        fid=filenum
    endif

    write (fid,'(/A32)') '================================'
    write (fid,'(A10,ES12.4,A10,/A10,ES12.4,A10,/A10,ES12.4,/A10,ES12.4)') 'F =', &
        this%F, 'V/nm', 'R =', this%R, 'nm', 'gamma =', this%gamma, 'theta_SC = ', &
        this%theta
    write (fid,'(A10,ES12.4,A10,/A10,ES12.4,A10)') 'W =', this%W, 'eV', &
            'kT =', this%kT, 'eV'
    write (fid,'(/A10,ES12.4,A10,/A10,ES12.4,A10,/A10,ES12.4,A10)') &
        'Jem =', this%Jem, 'A/nm^2', 'NotHeat =', this%heat,'W/nm^2', &
        '<DE> =', this%heat / this%Jem, 'eV'
    write (fid,'(/A10,ES12.4,A10,/A10,ES12.4,A10,/A10,ES12.4,A10,/A10,ES12.4)') &
        'xm =', this%xm, 'nm', 'xmax =', this%xmax, 'nm', &
        'Um =', this%Um, 'eV', 'Gamow =', this%Gam
    write (fid,'(/A10,ES12.4,A10,/A10,ES12.4)') 'L =', this%barlength, 'nm', &
            'ShapeFac =', this%Gam / (gg * sqrt(this%Um * this%barlength))
    write (fid,'(A10,ES12.4,A10,/A10,ES12.4,A10)') 'dG/dE@Ef =', this%maxbeta, &
            '(eV)^-1', 'dG/dE@Um =', this%minbeta, '(eV)^-1'
    write (fid,'(/A15,A20,/A15,A20,/A15,I20,/A15,A20,/A15,I20)') &
                                     'Regime:', trim(reg_str(this%regime)), &
                                     'Sharpness:', trim(sharp_str(this%sharpness)), &
                                     'Mode:', this%mode, &
                                     'Approximation:',trim(approx_str(this%approx)),&
                                     'Ierr:', this%ierr
    
    if (present(full) .and. full .and. allocated(this%xr)) then 
        !choose if V(x) is printed
        write(fid, '(/A32)') '  i        x_i          V(x_i)  '
        do i = 1, size(this%xr)
            write(fid, '(i3,ES14.4,ES15.5)') i, this%xr(i) , this%Vr(i)
        enddo
    endif
    
    write (fid,'(A32)') '--------------------------------'
end subroutine print_data


function fitpot(this)  result(var)
    !Fits V(x) data to F,R,gamma standard potential using L-V
    !minimization algorithm from SLATEC, implemented in the nlinfit subroutine here.
    
    type(EmissionData), intent(inout)   :: this
    
    real(dp)                :: p(3), F2, Fend, var, pmin(3), pmax(3)
    
    integer                 :: Nstart=3, i, fitinfo
    character(len=100)      :: Errormsg
    
    pmin = [0.d0, 1.d-5, 1.00001d0]
    pmax = [30.d0, 1.d4, 1.d5]
    
    p(1) = (this%Vr(Nstart) - this%Vr(Nstart-1)) &
            / (this%xr(Nstart) - this%xr(Nstart-1)) !guess for F (dV/dx @x=0)
    F2 = (this%Vr(Nstart+1) - this%Vr(Nstart)) / (this%xr(Nstart+1) -this%xr(Nstart))
    !F at the second point
    Fend = (this%Vr(size(this%Vr)) - this%Vr(size(this%Vr) - 1)) / &
            (this%xr(size(this%xr)) - this%xr(size(this%xr) - 1))!F at the end
    p(2) = abs(2.d0 / ((F2-p(1)) / (this%xr(Nstart) - this%xr(Nstart - 1))))
    !estimation for R
    p(3) = p(1) / Fend !estimation for gamma
    
    do i = 1, 3
        p(i) = min(p(i), pmax(i))
        p(i) = max(p(i), pmin(i))
    enddo
    
    if (debug > 2) print *, 'calling nlinfit. p =', p
    var = nlinfit(fun, this%xr, this%Vr, p, pmin, pmax, epsfit, fitinfo)
    
    if (fitinfo > 0 .and. fitinfo <=4 .and. var < varlim) then
        this%F = p(1)
        this%R = p(2)
        this%gamma = p(3)
    endif
    
    if (fitinfo == 0) call error_msg(this,'Wrong fitting input')
    if (fitinfo > 4) then
        write(Errormsg,*) 'Error in L-V fitting. INFO =', fitinfo
        call error_msg(this,Errormsg)
    endif

    contains

    pure function fun(x,p) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: p(:)
        real(dp)             :: y
        y = (p(1) * p(2) * x * (p(3) - 1.d0) + p(1) * x**2) &
            / (p(3) * x + p(2) * (p(3) - 1.d0))
    end function fun
end function fitpot


function fitpoly(this)  result(var)
    !Fits V(x) data to polynomial. Uses dpolfit from slatec.

    type(EmissionData), intent(inout)   :: this
    real(dp), dimension(size(this%xr))  :: ww, rr
    real(dp)                            :: var, eps
    integer                             :: ierr, N 
    
    ww(1) = -1.d0
    eps = epsfit !don't insert module parameter into slatec f77 function
    
    N = size(this%xr)
    
    if ((size(this%Apoly)) /=  (3*size(this%xr) + 3*Nmaxpoly + 3) .or. &
        .not.(allocated(this%Apoly))) then
        if (allocated(this%Apoly)) deallocate(this%Apoly)
        allocate (this%Apoly(3*size(this%xr) + 3*Nmaxpoly + 3))
    endif
    !make sure that the arrays are allocated properly and have the correct size
    if (debug > 1) print *, 'entering dpolfit. size(Apoly)=', size(this%Apoly)
    
    
    
    if (N<1 .or. Nmaxpoly > (N-1) .or. ww(1) /= -1.d0 ) then
        print *, 'Error in polynomial fit. N = ', N
        var = 1.d20
        return
    endif
    
     
    call dpolft(size(this%xr), this%xr, this%Vr, ww, Nmaxpoly, this%ndeg, eps,  &
                rr, ierr, this%Apoly)
    if(debug > 1) print *, 'exiting dpolfit' 
    var = eps
    
    if (debug > 1) print *, 'done polynomial fitting. ndeg=', this%ndeg, 'var =', var
    
    if (debug > 1 .and. ierr /= 1 .and. ierr /= 3) print *, &
        'error in polynomial fitting with ierr = ', ierr
            
end function fitpoly



!> Calculates the space charge reduction factor theta
function theta_SC(J, V, F) result(thetout)
    real(dp), parameter     :: k = 1.904e5
    real(dp), intent(in)    :: J, V, F !< Current density, total Voltage, local Field
    real(dp)                :: z, theta(2), thetout
    real(sp)                :: poly(3), work(8)
    complex(sp)             :: rt(2)
    integer                 :: ierr 
    
    z = k * sqrt(V) * J / F**2
    
    if (z < 1.e-3) then ! if z is too small the polynomial solution with sp fails
        thetout = 1.d0 - 4.d0 * z / 3.d0 + 3.d0 * z**2 !use asymptotic approximation
        !print *, 'zeta= ', z, 'theta= ', thetout
        return
    endif
    
    poly =  real([9. * z**2, -3.d0, -4. * z + 3.], sp)
    call RPQR79(2, poly, rt, ierr, work)
    
    if (ierr /= 0) then
        print *,'polynomial root finder in theta_SC returned error.'
        print *, 'poly = ', poly
    endif
    
    theta = real(rt,dp)
    
    if (abs(theta(1) - 2./3.) < abs(theta(2) - 2./3.)) then
        thetout = theta(1)
    else
        thetout = theta(2)
    endif
    !print *, 'zeta= ', z, 'theta= ', thetout

end function theta_SC




subroutine plot_barrier(this)
    !uses the pyplot module to plot the barrier defined by a specific Emisison object
    !this is useful mostly in debugging
    use pyplot_mod, only: pyplot
    use bspline, only: db1val, db1ink
    use std_mat, only: interp1
    
    integer, parameter      :: font = 32 , Nx = 512
    
    type (EmissionData), intent(inout)     :: this
    type (pyplot)           :: plt
    real(dp)                :: x(Nx), Ubar(Nx)
    real(dp)                :: Vplot(size(this%xr))
    integer                 :: ixrm, i, iflag, tempmode
    
    if (debug > 1) print *, "Entering plot_barrier"
    
    x = linspace(1.d-4,this%xmax,Nx)
    ixrm = size(this%xr)
    tempmode = this%mode
    
    
        
    if (this%mode == 1 .and. this%sharpness == BLUNT) then
        print *, 'changing mode'
        this%mode = 0
    endif
    
    do i =1, Nx
        Ubar(i) = bar(x(i))
        Ubar(i) = max(Ubar(i),-1.d0) 
    enddo
    
    call system("mkdir -p "//trim(outfolder))
    call system("mkdir -p "//trim(outfolder)//"/png "//trim(outfolder)//"/python")

    call plt%initialize(grid=.true.,xlabel='$x(nm)$',ylabel='$U(eV)$', &
            figsize=[20,15], font_size=font, title='FN-plot test', &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
    call plt%add_plot(x, Ubar, label='$barrier$', linestyle='b-', linewidth=2)
    if (allocated(this%xr)) then
        Vplot =this%W - this%Vr - Q / (this%xr + ( 0.5d0 * (this%xr**2)) / this%R)
        where (Vplot < -1.d0) Vplot = -1.d0
        call plt%add_plot(this%xr,Vplot,label='$barrier points$', &
                    linestyle='r*',markersize = 10)
    endif
    
    call plt%savefig(trim(outfolder)//'/png/barrierplot.png', &
            pyfile=trim(outfolder)//'/python/barrierplot.plot.py')
    
    this%mode = tempmode
    

    contains
    
    function bar(x) result(V)!sphere barrier model
        real(dp), intent(in)    :: x
        real(dp)                :: V, Vinterp(2)
        integer                 :: iflag, inbvx
        
        select case (this%mode)
            case (2) !interpolation mode.. Use splines
                inbvx = 1
                call db1val(x, idx, this%tx, size(this%Vr), knotx, &
                        this%bcoef, Vinterp(1), iflag, inbvx)
                if (iflag /= 0) then !out of bounds.. do linear extrapolation 
                    Vinterp(1) = this%Vr(ixrm) + (x - this%xr(ixrm)) * &
                    (this%Vr(ixrm)-this%Vr(ixrm-1)) / (this%xr(ixrm)-this%xr(ixrm-1))
                endif
                V = this%W - Vinterp(1) - Q / (x + ( 0.5d0 * (x**2)) / this%R)
            case (-22) ! fitted polynomial mode. Evaluate fitted polynomial
                if (x<= this%xr(ixrm)) then
                    call dp1vlu(this%Ndeg, 0, x, Vinterp(1), Vinterp(2:), this%Apoly)
                else !out of bounds. Polynomial misbehaves. Use extrapolation
                    call dp1vlu(this%Ndeg, 1, this%xr(ixrm), & !polyval + der in xrm
                            Vinterp(1), Vinterp(2:), this%Apoly)
                    Vinterp(1) = Vinterp(1) + Vinterp(2) * (x - this%xr(ixrm))
                    !linear extrapolation using the derivative of the polynomial
                    Vinterp(1) = Vinterp(1)
                endif
                V = this%W - Vinterp(1) - Q / (x + ( 0.5d0 * (x**2)) / this%R)
            case (0,-12, -1)
            !use standard F,R,gamma model for the electrostatic potential
                V = this%W - (this%F * this%R * x*(this%gamma - 1.d0)  &
                + this%F * x**2) / (this%gamma * x + this%R * (this%gamma - 1.d0)) &
                - Q / (x + ( 0.5d0 * (x**2)) / this%R)
            case default
                print *, 'Error: invalid mode in bar(). Mode = ', this%mode  
                stop
        end select
    end function bar
end subroutine plot_barrier


subroutine desetroy(this)
!frees all the memory allocated for an Emission type "object"
    type(EmissionData), intent(inout)   :: this
    
    if (this%mode == 2) then
        this%mode = 1 !return to initial mode 1
        deallocate(this%bcoef, this%tx)
    endif
    
    if (allocated(this%xr)) deallocate(this%xr,this%Vr)
    if (allocated(this%Apoly)) deallocate(this%Apoly)
    
end subroutine desetroy


subroutine C_wrapper(passdata, ifun) bind(c)
! a wrapper to be able to call getelec functions from C language
! passdata is a C struct corresponding to the EmissionData type
! ifun gives which getelec subroutine will be called on the specific data struct
! case ifun = 0: cur_dens(), 
! case ifun = 1: print_data(),
! case ifun = 2: print_data(full = .true.)
! case ifun = 3: plot_barrier() 
! case ifun = 4: theta_SC()
! case ifun = 5: cur_dens_SC()   

    use iso_c_binding  
                                                  
    type, bind(c)   :: cstruct
        real(c_double)      :: F, W, R, gamma, Temp !input parameters
        real(c_double)      :: Jem, heat !ouput parameters
        type(c_ptr)         :: xr, Vr     !input vectors
        integer(c_int)      :: regime, sharp  !output chars showing regimes
        integer(c_int)      :: Nr, approx, mode, ierr !len of vectors xr, Vr and full
        real(c_double)      :: voltage, theta
    end type cstruct

    type(cstruct), intent(inout)    :: passdata
    integer(c_int), intent(in), value :: ifun
   
    real(c_double), pointer         :: xr_fptr(:), Vr_fptr(:)
    type(EmissionData)              :: this
    integer                         :: i
    
    if (ifun == 4) then
        passdata%Temp = theta_SC(passdata%F, passdata%voltage, passdata%R)
        return
    endif

    
    !copy the members of the c struct to the fortran type
    this%F = passdata%F; this%W = passdata%W; this%R = passdata%R; 
    this%kT = kBoltz * passdata%Temp; this%gamma = passdata%gamma; 
    this%mode = passdata%mode; this%approx = passdata%approx
    this%Jem = passdata%Jem; this%heat = passdata%heat
    this%ierr = passdata%ierr; this%regime = passdata%regime
    this%sharpness = passdata%sharp
    
    
    if (this%mode /= 0 .and. this%mode /= -1) then
        if (passdata%Nr == 0) then
            print *, 'Nr = ', passdata%Nr
            call print_data(this)
            stop 'Error: C_wrapper. Incompatible mode with length of potential array'
        else
            call c_f_pointer(passdata%xr, xr_fptr, [passdata%Nr])
            !copy pointers to fortran from c
            call c_f_pointer(passdata%Vr, Vr_fptr, [passdata%Nr])
            
            allocate(this%xr(passdata%Nr), this%Vr(passdata%Nr))
            
            do i = 1, passdata%Nr  !copy c input data to object data
                this%xr(i) = xr_fptr(i)
                this%Vr(i) = Vr_fptr(i)
            enddo
        endif
    endif

    if (ifun >= 5) then
        call cur_dens_SC(this, passdata%voltage)
    else
        call cur_dens(this)
    endif

    passdata%Jem = this%Jem
    passdata%heat = this%heat
    passdata%regime = this%regime
    passdata%sharp = this%sharpness
    passdata%ierr = this%ierr
    passdata%theta = this%theta

    select case (ifun)
        case (0)
            continue
        case (1)
            call print_data(this)
        case (2)
            call print_data(this,.true.)
        case (3)
            call plot_barrier(this)
        case default
            return
    end select

    call desetroy(this)

end subroutine C_wrapper


function nlinfit(fun, xdata, ydata, p0, pmin, pmax, tol, info) result(var)
 
!Fit data to arbitrary function by using the Levenberg-Marquardt algorithm
!implemented in the dnls1e function in slatec. This implementation gives the option
!to give extreme values for the fitted parameters and introduce a y-shift to match
!the first value of the data. Handy in case there is an arbitrary adjustable additive
!parameter as in the case of FN-plot analysis with the sigma prefactor.

    real(dp), intent(in)        :: xdata(:), ydata(:), pmin(:), pmax(:), tol 
    !input x, y data, vectors with min and max params, tolerance
    real(dp), intent(inout)     :: p0(:)    !in: initial guess, out:result
    integer, intent(out)        :: info
    
    interface                               !user provided function to be fit
        function fun(x,p) result(y)
            implicit none
            integer,parameter    :: dp = 8
            real(8), intent(in) :: x
            real(8), intent(in) :: p(:)
            real(dp)             :: y
        end function fun
    end interface

    real(dp)                    :: var, fvec(size(xdata)), yshiftloc
    integer                     :: i, m, n, iwa(size(p0)), lwa, iopt, nprint
    integer, save               :: Nvals = 0
    real(dp)                    :: wa(size(p0) * (size(xdata) + 5) + size(xdata))
    
    m = size(xdata)
    n = size(p0)
    lwa = n *(m + 5) + m
    iopt = 1
    nprint = 0
    
    if (debug > 2) print *, 'Calling dnlse1'

    call DNLS1E(fcn, iopt, m, n, p0, fvec, tol, nprint, info, iwa, wa, lwa)
    if (info > 3 .or. info == 0) then
        var = 1.d100
    else
        var = sum(fvec) / sum(abs(ydata)) !relative error (mean error/mean |values|)
    endif

    if (debug > 1) print *, 'Model fitting completed. info =', &
                info, 'tol = ', tol, 'Nvals = ', Nvals, 'var =', var
    
    contains

    subroutine fcn(iflag, m, n, p, fvec, fjac, ldfjac)
        
        integer, intent(in)     :: m, n
        real(dp), intent(in)    :: p(n)
        real(dp), intent(inout) :: fvec(m), fjac
        integer, intent(inout)  :: iflag, ldfjac
        
        real(dp)                :: peval(n)
        real(dp)                :: multiplier, funeval
        
        do i = 1, n ! limit peval between pmin and pmax
            peval(i) = max(p(i), pmin(i))
            peval(i) = min(peval(i), pmax(i))
        enddo
        
        multiplier  = sum(abs(peval - p))
        do i = 1, m
            funeval = fun(xdata(i), peval)
            fvec(i) = abs(funeval - ydata(i))
        enddo

        if (multiplier > 1.d-10) &
            fvec = fvec * 1.d2 * exp(multiplier)
            
        Nvals = Nvals + 1
        
    end subroutine fcn

end function nlinfit


subroutine read_params()

    character(len=13)   :: str
    logical             :: ex
    integer             :: fid = fidparams
    character(len=40), parameter  :: error = &
        'ERROR: GetelecPar.in is not correct.'
    
    inquire(file=paramfile, exist=ex)
    if (.not. ex) then
        print *, 'GETELEC: Parameters input file not found. Default values used.'
        readparams = .true.
    endif
    
    if (readparams) return
    
    open(fid, file = paramfile, action = 'read')
    read(fid,'(2A)')
    
    read(fid,*) str, debug
    if (str(1:5)/='debug') stop error 
    if (debug > 1)  print *, 'debug read', debug

    read(fid,*) str, xlim
    if (str(1:4)/='xlim') stop error
    if (debug > 1)  print *, 'xlim read', xlim
    
    read(fid,*) str, nlimfield
    if (str(1:9)/='nlimfield') stop error
    if (debug > 1)  print *, 'nlimfield read', nlimfield
    
    read(fid,*) str, nlimthermal
    if (str(1:11)/='nlimthermal') stop error
    if (debug > 1)  print *, 'nlimthermal read', nlimthermal
    
    read(fid,*) str, nmaxlim
    if (str(1:7)/='nmaxlim') stop error
    if (debug > 1)  print *, 'nmaxlim read', nmaxlim
    
    read(fid,*) str, gammalim
    if (str(1:8)/='gammalim') stop error
    if (debug > 1)  print *, 'gammalim read', gammalim
    
    read(fid,*) str, varlim
    if (str(1:6)/='varlim') stop error
    if (debug > 1)  print *, 'varlim read', varlim
    
    read(fid,*) str, epsfit
    if (str(1:6)/='epsfit') stop error
    if (debug > 1)  print *, 'epsfit read', epsfit
    
    read(fid,*) str, Nmaxpoly
    if (str(1:8)/='Nmaxpoly') stop error
    if (debug > 1)  print *, 'Nmaxpoly read', Nmaxpoly
    
    read(fid,*) str, spectra
    if (str(1:7)/='spectra') stop error 
    if (debug > 1)  print *, 'spectra read', spectra
    
    read(fid,*) str, above
    if (str(1:7)/='aboveUm') stop error 
    if (debug > 1)  print *, 'aboveUm read', above
    
    read(fid,*) str, outfolder
    if (str(1:6)/='output') stop error
    outfolder = trim(outfolder) 
    if (debug > 1)  print *, 'outfolder read', outfolder
    
    close(fid)
    
    if (debug > 2)  print *, "Parameters read"
    
    readparams = .true.

end subroutine read_params


subroutine error_msg(this, msg)
    ! prints error message to the error file
    type(EmissionData), intent(inout)  :: this
    character(len=*)          :: msg
    
    if (debug == 0) return
    
    open(fiderr, file = trim(outfolder)//'/'//errorfile, action ='write', &
                    access ='append')
    call print_data(this, .true., fiderr)
    write(fiderr,*) trim(msg)
    close(fiderr)
end subroutine error_msg

    
end module GeTElEC
