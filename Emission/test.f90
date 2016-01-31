program test

use emission, only: elec_emit

double precision::J

J=elec_emit(5.d0,4.5d0,5.d0,0.d0)

print *, J

end program