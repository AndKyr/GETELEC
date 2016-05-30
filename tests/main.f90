program main

use emission, only: gamow_general, EmissionData, print_data, J_num_integ, cur_dens, kBoltz
use pyplot_mod, only: pyplot
use std_mat, only: linspace

implicit none

integer,parameter       :: dp=8, Nf=1024, font=35
real(dp), parameter     :: Fmin=0.3d0, Fmax=7.d0

real(dp)                :: T=800.d0, Fi(Nf), Ji(Nf), heati(Nf), t1,t2, Japp(Nf), heatapp(Nf)
!real(dp), allocatable   :: Jib(Nf) , Jis(Nf)
integer                 :: i , fidout = 1564
character               :: reg(Nf), sharp(Nf), regap(Nf), sharpap(Nf)

type(EmissionData)      :: new
type(pyplot)            :: plt 


new%kT=kBoltz*T
new%R = 2.d0
new%gamma = 10.d0
Fi=linspace(Fmin,Fmax,Nf)
call cpu_time(t1)
do i=1,Nf
    
    new%F=Fi(i)
    
    new%full = .true.
    call cur_dens(new)
    Ji(i) = new%Jem
    heati(i) = new%heat
    reg(i) = new%regime
    sharp(i) = new%sharpness

    new%full = .false.
    call cur_dens(new)
    Japp(i) = new%Jem
    heatapp(i) = new%heat
    regap(i) = new%regime
    sharpap(i) = new%sharpness
    
enddo
call cpu_time(t2)

open(fidout,file="python/main.py",action="write",status="replace")
write(fidout,*)'J = np.array([', (Ji(i), ',', i=1,Nf-1), Ji(Nf),'])'
write(fidout,*)'F = np.array([', (Fi(i), ',', i=1,Nf-1), Fi(Nf),'])'
write(fidout,*)'Japp = np.array([', (Japp(i), ',', i=1,Nf-1), Japp(Nf),'])'
write(fidout,*)'heat = np.array([', (heati(i), ',', i=1,Nf-1), heati(Nf),'])'
write(fidout,*)'heatapp = np.array([', (heatapp(i), ',', i=1,Nf-1), heatapp(Nf),'])'
write(fidout,*)'reg = np.array([', (ichar(reg(i)), ',', i=1,Nf-1), ichar(reg(Nf)),'])'
write(fidout,*)'sharp = np.array([', (ichar(sharp(i)), ',', i=1,Nf-1), ichar(sharp(Nf)),'])'
write(fidout,*)'regap = np.array([', (ichar(regap(i)), ',', i=1,Nf-1), ichar(regap(Nf)),'])'
write(fidout,*)'sharpap = np.array([', (ichar(sharpap(i)), ',', i=1,Nf-1), ichar(sharpap(Nf)),'])'

close(fidout)

call plt%initialize(grid=.true.,xlabel='$1/F [nm/V]$',ylabel='$J (A/nm^2)$', &
            figsize=[20,15], font_size=font, title='FN-plot test', &
            legend=.true.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
call plt%add_plot(1/Fi,log10(Ji),label='$current$', &
                    linestyle='b-',linewidth=2)
                    
call plt%add_plot(1/Fi,log10(abs(heati)),label='$heat$', &
                    linestyle='r-',linewidth=2)
                    
call plt%add_plot(1/Fi,log10(Japp),label='$current-approx$', &
                    linestyle='b--',linewidth=2)
                    
call plt%add_plot(1/Fi,log10(abs(heatapp)),label='$heat-approx$', &
                    linestyle='r--',linewidth=2)
                    
call plt%savefig('png/FNplot.png', pyfile='python/FNplot.py')

call print_data(new)
end program
