program main

use GeTElEC, only: EmissionData, cur_dens, kBoltz
use pyplot_mod, only: pyplot
use std_mat, only: linspace

implicit none

integer,parameter       :: dp=8, Nf=256, font=40, lw = 3
real(dp), parameter     :: Fmin=3.d0, Fmax=10.d0

real(dp), dimension(Nf) :: Fi, Umi, betamin, betamax, error
integer                 :: i, j

real(dp)                :: Jfull, Jgtf, Wi(3) = [4.2d0, 4.6d0, 5.2d0]


type(EmissionData)      :: this
type(pyplot)            :: plt

call plt%initialize(grid=.true.,xlabel='$error$', &
            figsize=[20,15], font_size=font, &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
        
this%gamma = 10.d0

this%approx = 0
this%R = 1.d2
this%mode = 0

Fi=linspace(Fmin, Fmax, Nf)

do j = 1,3
    this%W = Wi(j)
    do i = 1, Nf
        this%approx = 1
        this%F= Fi(i)
        call cur_dens(this)
        Jfull = this%Jem
        Umi(i) = this%Um
        betamax(i) = this%maxbeta
        betamin(i) = this%minbeta
        
        this%approx = 0
        call cur_dens(this)
        Jgtf = this%Jem
        
        error(i) = 100. * abs(Jfull - Jgtf) / Jfull
    enddo
        
    call plt%add_plot(error, 1./(Umi * betamax), linestyle='b-',label='Um',linewidth=2)
!    call plt%add_plot(error, 1/log(betamax / betamin),linestyle='r-',label='dbeta',linewidth=2)
    !call plt%add_plot(Fi,betamax-betamin,linestyle='k-',label='betadiff',linewidth=2)
enddo

call plt%savefig('png/highFplot.png', pyfile='python/highFplot.plot.py')                    

end program
