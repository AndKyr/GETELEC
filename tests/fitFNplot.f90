program fitFN

use GeTElEc, only: fitFNplot, dp, EmissionData, cur_dens, print_data, kBoltz
use pyplot_mod, only: pyplot

integer, parameter      :: fid = 549687, Nmaxeval = 10000, font = 30

integer                 :: Ndata

real(dp), allocatable   :: xdata(:), ydata(:), ycalc(:)
real(dp)                :: var, epsfit = 1.d-12, params(5), pmin(5), pmax(5), yshift
type(pyplot)            :: plt



open(fid, file = 'IVplot.dat', action = 'read')

read (fid, *) Ndata
read (fid, '(A)')

allocate(xdata(Ndata), ydata(Ndata), ycalc(Ndata))

do i = 1, Ndata
    read (fid, *) xdata(i), ydata(i)
enddo

pmin =   [0.5d0 * sum(xdata)/Ndata, 4.5d0, 0.1d0, 10.d0, 290.d0]
params = [3 .d0 * sum(xdata)/Ndata, 4.5d0, 2.d0, 10.d0, 300.d0]
pmax =   [20.d0 * sum(xdata)/Ndata, 4.5d0, 1.d2, 10.d0, 310.d0]


var = fitFNplot(xdata, ydata, params, pmin, pmax, epsfit, Nmaxeval, yshift)

print *, 'result: [beta, W, R, gamma, T]: ', params
print *, 'variance:', var, 'yshift = ', yshift
print *, 'fbounds:', params(2)/minval(xdata), params(2) / maxval(xdata) 


call plt%initialize(grid=.true.,xlabel='$1/F \quad [nm/V]$',ylabel='$log(I) [A]$', &
        figsize=[20,15], font_size=font, title='FN-plot test', &
        legend=.true.,axis_equal=.false., legend_fontsize=font, &
        xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
call plt%add_plot(xdata, ydata, label='$data$', linestyle='b*', linewidth=2)

do i = 1, Ndata
    ycalc(i) = fun(xdata(i), params)
enddo

call plt%add_plot(xdata, ycalc + yshift, label='$fitted$', linestyle='r-', linewidth=2)

call plt%savefig('barrierplot.png', pyfile='barrierplot.plt.py')


    contains

    function fun(x, p) result(y)
    
        real(dp), intent(in)    :: x, p(:)
        real(dp)                :: y, dpmin, dpmax
        type(EmissionData)          :: this
         
        
        this%F = p(1) / x !F = beta * V  = beta / xdata
        this%W = p(2) !work function
        this%R = p(3) !radius
        this%gamma = p(4)  ! gamma
        this%kT = kBoltz * p(5) !temperature
        
        call cur_dens(this)
        y = log(this%Jem)
    end function fun

end program




