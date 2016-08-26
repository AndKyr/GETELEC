program KXerror

use GeTElEc, only: dp, EmissionData, gamow_num, gamow_KX, kBoltz
use std_mat, only: linspace
use pyplot_mod, only: pyplot

integer, parameter      :: fid = 549687, Nmaxeval = 10000, font = 38, Nchi = 256
real(dp), parameter     :: chimin = 0.09, chimax = 0.11, gammai(2) = [10.d0, 20.d0]
real(dp), parameter     :: Wconst = 4.5d0, Fconst = 3.d0, Rconst = 10.d0

real(dp),dimension(Nchi):: chi, Wi, Fi, Ri, GKX, Gnum, Gerr, DKX, Dnum   
type(pyplot)            :: plt
type(EmissionData)      :: thisKX, thisnum
integer                 :: i, j, info

character(len=3)        :: lstyle(4)


chi = linspace(chimin, chimax, Nchi)
lstyle = ['b- ','r--','m-.','k: ']

call plt%initialize(grid=.true.,xlabel='$\chi$',ylabel='$\Delta G / G (\%)$', &
        figsize=[20,15], font_size=font, title='G error', &
        legend=.true.,axis_equal=.false., legend_fontsize=font, &
        xtick_labelsize=font,ytick_labelsize=font,axes_labelsize=font)
        
do j=1, 4
    select case(j)
        case (1)
            Wi = Wconst
            Ri = Rconst
            Fi = Wi / (chi * Ri)
            gamma = gammai(1)
        case (2)
            Wi = Wconst
            Fi = Rconst
            Ri = Wi / (chi * Fi)
            gamma = gammai(1)
        case (3)
            Ri = Rconst
            Fi = Fconst
            Wi = chi * Fi * Ri
            gamma = gammai(1)
        case (4)
            Wi = Wconst
            Fi = Rconst
            Ri = Wi / (chi * Fi)
            gamma = gammai(2)
    end select
    
    do i = 1, Nchi
        thisnum%F = Fi(i)
        thisnum%W = Wi(i)
        thisnum%R = Ri(i)
        thisnum%gamma = gamma
        thisnum%mode = 0
        thisKX = thisnum
        info = gamow_KX(thisKX, .false.)
        call gamow_num(thisnum, .true.)
        GKX(i) = thisKX%Gam
        Gnum(i) = thisnum%Gam
    enddo
    
    DKX = 1.d0 / (1.d0 + exp(GKX))
    Dnum = 1.d0 / (1.d0 + exp(Gnum))
    
    Gerr = 100.d0* abs(Gnum - GKX) / Gnum
            
    call plt%add_plot(chi, Gerr, label='$data$', linestyle=lstyle(j), linewidth=2)

enddo

call plt%savefig('png/Gerror.png', pyfile='python/Gerror.py')

end program



