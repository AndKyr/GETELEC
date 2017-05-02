program main

use GeTElEC, only: gamow_general, EmissionData, print_data, cur_dens, kBoltz
use pyplot_mod, only: pyplot
use std_mat, only: linspace

implicit none

integer,parameter       :: dp=8, Nf=512, font=40, lw = 3
real(dp), parameter     :: Fmin=0.4d0, Fmax=14.d0, T = 2000.0

real(dp), dimension(Nf) :: Fi, Jfs, Jis, Jts, Jfb, Jib, Jtb, &
                           hfs, his, hts, hfb, hib, htb, Japp, happ, &
                           Ffs, Fis, Fts, Ffb, Fib, Ftb
integer                 :: i, j, fs=0, is=0, ts=0, fb=0, ib=0, tb=0
real(dp)                :: Ri(3)

logical, parameter      :: DE = .true.

type(EmissionData)      :: this
type(pyplot)            :: plt1, plt2, plt3

call plt1%initialize(grid=.true.,xlabel='$1/F \ [nm/V]$',ylabel='$J \ [A/nm^2]$', &
            figsize=[20,15], font_size=font, &
            legend=.false.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
call plt2%initialize(grid=.true.,xlabel='$1/F \ [nm/V]$',ylabel='$P_N \ [W/nm^2]$', &
            figsize=[20,15], font_size=font, &
            legend=.false.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
            
if (DE) call plt3%initialize(grid=.true.,xlabel='$1/F \ [nm/V]$', &
            ylabel='$\Delta E [eV]$', figsize=[20,15], font_size=font, &
            legend=.false.,axis_equal=.false., legend_fontsize=font, &
            xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)


this%kT=kBoltz*T
this%gamma = 15.d0
!this%W = 4.5d0
Fi=1.d0/linspace(1.d0/Fmax,1.d0/Fmin,Nf)
Ri = [1.d0, 5.d0, 200.d0]
this%mode = 0




do j = 1,3
        fs=0
        is=0
        ts=0
        fb=0
        ib=0
        tb=0
        this%R = Ri(j)
    do i=1,Nf

        
        this%F= Fi(i)
        this%full = .true.
        this%mode = 0
        call cur_dens(this)
        
        if (this%regime == 'F') then
            if (this%sharpness == 'S') then
                fs = fs + 1
                Jfs(fs) = this%Jem
                hfs(fs) = this%heat
                Ffs(fs) = this%F
            else
                fb = fb + 1
                Jfb(fb) = this%Jem
                hfb(fb) = this%heat
                Ffb(fb) = this%F
            endif
        elseif (this%regime == 'I') then
            if (this%sharpness == 'S') then
                is = is + 1
                Jis(is) = this%Jem
                his(is) = this%heat
                Fis(is) = this%F
            else
                ib = ib + 1
                Jib(ib) = this%Jem
                hib(ib) = this%heat
                Fib(ib) = this%F
            endif
        else
            if (this%sharpness == 'S') then
                ts = ts + 1
                Jts(ts) = this%Jem
                hts(ts) = this%heat
                Fts(ts) = this%F
            else
                tb = tb + 1
                Jtb(tb) = this%Jem
                htb(tb) = this%heat
                Ftb(tb) = this%F
            endif
        endif
    if (this%R > 0.d0) then
        this%full = .false.
        call cur_dens(this)
        Japp(i) = this%Jem
        happ(i) = this%heat
    endif
        
    enddo
    
    if (this%R < 1.d2) then
              
        if (fs > 0) call plt1%add_plot(1./Ffs(1:fs), Jfs(1:fs), linestyle='b-', &
                        label='fs', linewidth=lw, yscale = 'log')
        if (fb > 0) call plt1%add_plot(1./Ffb(1:fb), Jfb(1:fb), linestyle='b--', &
                        yscale = 'log', label='fb', linewidth=lw)
        if (ib > 0) call plt1%add_plot(1./Fib(1:ib), Jib(1:ib), linestyle='k--', &
                        yscale = 'log', label='ib', linewidth=lw)                    
        if (is > 0) call plt1%add_plot(1./Fis(1:is), Jis(1:is), linestyle='k-', &
                        yscale = 'log', label='is', linewidth=lw)
        if (tb > 0) call plt1%add_plot(1./Ftb(1:tb), Jtb(1:tb), linestyle='r--', &
                        yscale = 'log', label='tb', linewidth=lw)
        if (ts > 0) call plt1%add_plot(1./Fts(1:ts), Jts(1:ts), linestyle='r-', &
                        yscale = 'log', label='ts', linewidth=lw)

                              
        if (fs > 0) call plt2%add_plot(1./Ffs(1:fs), abs(hfs(1:fs)), linestyle='b-',&
                        label='fs', linewidth=lw, yscale = 'log')
        if (fb > 0) call plt2%add_plot(1./Ffb(1:fb), abs(hfb(1:fb)),linestyle='b--',&
                        yscale = 'log', label='fb', linewidth=lw)
        if (ib > 0) call plt2%add_plot(1./Fib(1:ib), abs(hib(1:ib)),linestyle='k--',&
                        yscale = 'log', label='ib', linewidth=lw)                    
        if (is > 0) call plt2%add_plot(1./Fis(1:is), abs(his(1:is)), linestyle='k-',&
                        yscale = 'log', label='is', linewidth=lw)
        if (tb > 0) call plt2%add_plot(1./Ftb(1:tb), abs(htb(1:tb)),linestyle='r--',&
                        yscale = 'log', label='tb', linewidth=lw)
        if (ts > 0) call plt2%add_plot(1./Fts(1:ts), abs(hts(1:ts)), linestyle='r-',&
                        yscale = 'log', label='ts', linewidth=lw)
    endif
                    
    if (this%R > 0.d0) then
        call plt1%add_plot(1/Fi,Japp,linestyle='m:', &
                        label='app',linewidth=2, yscale = 'log')
        call plt2%add_plot(1./Fi, abs(happ), linestyle='m:', label='app', &
                    linewidth=2, yscale = 'log')
    endif
    
    if (DE .and. this%R > 100.d0) call plt3%add_plot(1/Fi, happ/Japp, &
                                    linestyle='k-', label='app', linewidth=2)
        
        
enddo

call plt1%savefig('png/Jplot.png', pyfile='python/Jplot2.plot.py')                    
call plt2%savefig('png/heatplot.png', pyfile='python/heatplot2.plot.py')
if (DE) call plt3%savefig('png/DE.png', pyfile='python/DEplot.plot.py')

end program
