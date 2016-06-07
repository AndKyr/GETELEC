program test

use std_mat

real(dp)    :: x , xi(32), yi(32), y(2)

xi = linspace(0.d0,1.d0,32)
yi =  xi**2 + 3.d0 * xi

x = 0.1598752
y = interp1(xi,yi,xi(8))
print *, yi(8), y
!print *,   x**2 + 3.d0 * x

end program










	
