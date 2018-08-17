program current

use GeTElEC, only: EmissionData, dp, cur_dens, print_data, plot_barrier, debug, cur_dens_SC
use std_mat, only: linspace

type(EmissionData)      :: this
integer                 :: i, Nx = 32
character(len=32)       :: arg
real(dp)                :: t1,t2, kBoltz = 8.6173324d-5, voltage
    
i = iargc()      
if (i < 3) then
    print *, 'Give at least 3 floats. Field in V/nm, work function  in eV &
    and temperature in degrees Kelvin'
    print *, 'Optionally give afterwards R, mode, approximation type &
    and gamma'
    print *, "Approximation types are: 1: Full, 0: GTF, -1: FN, -2: RLD "
    print *, "If you want space charge included give voltage [V]"
    stop
endif


call getarg(1, arg)
arg= trim(arg)
read(arg,*) this%F

call getarg(2, arg)
arg= trim(arg)
read(arg,*) this%W

call getarg(3, arg)
arg= trim(arg)
read(arg,*) this%kT
this%kT = this%kT*kBoltz

if (i >= 4) then
    call getarg(4, arg)
    arg= trim(arg)
    read(arg,*) this%R
endif

if (i >= 5) then
    call getarg(5, arg)
    arg= trim(arg)
    read(arg,*) this%mode
endif

if (i >= 6) then
    call getarg(6, arg)
    read(arg,*) this%approx
endif
if (i >= 7) then
    call getarg(7, arg)
    arg= trim(arg)
    read(arg,*) this%gamma
endif
if (i >= 8) then
    call getarg(8, arg)
    arg= trim(arg)
    read(arg,*) voltage
endif


if (this%mode > 0 .or. this%mode < -1) then
    allocate(this%xr(Nx), this%Vr(Nx))
    this%xr = linspace(0.d0,3.d0,Nx)
    this%Vr = (this%F * this%R * this%xr*(this%gamma - 1.d0) + this%F * this%xr**2) &
                / (this%gamma * this%xr + this%R * (this%gamma - 1.d0))
endif
call cpu_time(t1)
if (i>=8) call cur_dens_SC(this, voltage)
if (i<8) call cur_dens(this)
call cpu_time(t2)
call print_data(this,.true.)
if (debug > 0) then
    print *, 'Timings: ', this%timings
    print *, 'total elapsed time:', t2-t1
endif


end program current



