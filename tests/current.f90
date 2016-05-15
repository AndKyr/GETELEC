program current

use emission, only: EmissionData, dp, cur_dens, kBoltz, print_data

type(EmissionData)      :: this
integer                 :: i
character(len=32)       :: arg
    
i = iargc()      
if (i < 3) then
    print *, 'Give at least 3 floats. Field, work function and temperature'
    stop
endif


call getarg(1, arg)
arg= trim(arg)
print *, arg
read(arg,*) this%F

call getarg(2, arg)
arg= trim(arg)
read(arg,*) this%W

call getarg(3, arg)
arg= trim(arg)
read(arg,*) this%kT
this%kT = this%kT*kBoltz

if (i == 4) then
    call getarg(4, arg)
    arg= trim(arg)
    read(arg,*) this%R
endif
if (i == 5) then
    call getarg(5, arg)
    arg= trim(arg)
    read(arg,*) this%gamma
endif

call cur_dens(this)
call print_data(this)

end program current



