program test


use emission
implicit none

type(EmissionData)      :: this

this%kT=0.07d0
this%full = .false.
this%F = 2.d0

call cur_dens(this)

call print_data(this)

this%full = .true.

call cur_dens(this)

call print_data(this) 



end program
