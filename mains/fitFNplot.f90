program fitFN

use GeTElEc, only fitFNplot

integer, parameter      :: fid = 549687

integer                 :: Ndata

real(dp), allocatable   :: xdata(:), ydata(:)
real(dp)                :: var, betas = 1.d0, works = 4.5d0, radii = 5.d0, &
                            gammas = 10.d0, Temps = 400.d0



open(fid, file = 'IVplot.dat', action = 'read')

read (fid, *) Ndata
read (fid, '(A)')

allocate(xdata(Ndata), ydata(Ndata))

do i = 1, Ndata
    read (fid, *) xdata(i), ydata(i)
enddo


var = fitFNplot(xdata,ydata,betas, works, radii, gammas, Temps)

print *, 'result: [beta, W, R, gamma, T]: ', betas, works, radii, gammas, Temps
print *, 'variance:', var

end program




