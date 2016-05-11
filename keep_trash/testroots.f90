program testroot

integer             :: ndeg=2, ierr
real                :: p(3), work(8)
complex             :: roots(2)

p = [1., -2., 4.]
call RPQR79(ndeg,p,roots,ierr,work)

print *, 'Rt=', real(roots), '|  ierr = ', ierr

end program testroot