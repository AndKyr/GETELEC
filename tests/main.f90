program main

use emission, only: gamow_general, EmissionData, print_data, J_num_integ, cur_dens, kBoltz
use pyplot_mod, only: pyplot
use std_mat, only: linspace

implicit none

integer,parameter   :: dp=8, Nf=256, font=35
real(dp), parameter :: Fmin=0.05d0, Fmax=7.d0

real(dp)            :: T=700.d0, Fi(Nf), Ji(Nf), heati(Nf), t1,t2
integer             :: i

type(EmissionData)      :: old,new
type(pyplot)            :: plt 



old%kT=kBoltz*T
old%R = 10.d0
Fi=linspace(Fmin,Fmax,Nf)
new=old
call cpu_time(t1)
do i=1,Nf
    
    new%F=Fi(i)
    
    call cur_dens(new)
    Ji(i) = new%Jem
    heati(i) = new%heat
    print *, Ji(i), heati(i), new%regime, new%sharpness
enddo
call cpu_time(t2)
print *, 'Elapsed time:', t2-t1

call plt%initialize(grid=.true.,xlabel='$1/F [nm/V]$',ylabel='$J (A/nm^2)$', &
            figsize=[20,15], font_size=font, title='FN-plot test', &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
call plt%add_plot(1/Fi,log10(Ji),label='$current$', &
                    linestyle='b-',linewidth=2)
                    
call plt%add_plot(1/Fi,log10(abs(heati)),label='$heat$', &
                    linestyle='r-',linewidth=2)
                    
call plt%savefig('FNplot.png', pyfile='FNplot.py')

call print_data(new)
end program
