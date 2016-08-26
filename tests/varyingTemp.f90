program main

use GeTElEC, only: EmissionData, cur_dens, kBoltz
use pyplot_mod, only: pyplot
use std_mat, only: linspace

implicit none

integer,parameter       :: dp=8, Nf=2048, font=40, lw = 3
real(dp), parameter     :: Tmin=4.d0, Tmax=4.d3

real(dp), dimension(Nf) :: Ti, DEi
integer                 :: i, j
real(dp)                :: R, Fi(3)

logical, parameter      :: DE = .true.
character(len=2)        :: ls(3)

type(EmissionData)      :: this
type(pyplot)            :: plt

call plt%initialize(grid=.true.,xlabel='$T \ [K]$',ylabel='$\Delta E \ [eV]$', &
            figsize=[20,15], font_size=font, &
            legend=.false.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
        
this%gamma = 15.d0


Fi = [4.d0, 6.d0, 8.d0]
ls = ['b-','g-','r-']
this%full = .false.
this%R = 1.d2
this%mode = 0

Ti=linspace(Tmin, Tmax, Nf)
do j = 1,3
    this%F= Fi(j)
    do i=1,Nf

        this%kT = kBoltz * Ti(i)
        
        call cur_dens(this)
        DEi(i) = -this%heat / this%Jem
    enddo
        
    call plt%add_plot(Ti,DEi,linestyle=ls(j),label='F=',linewidth=2)
enddo
    

call plt%savefig('png/varytemp.png', pyfile='python/varytemp.plt.py')                    

end program
