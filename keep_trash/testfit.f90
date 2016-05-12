program testfit  !Fits V(x) data to F,R,gamma standard potential using L-V
    !minimization module
    use pyplot_mod, only: pyplot
    use fitting, only: sp, mrqmin, i4b

    integer(I4B), parameter      :: ndata = 32, Np = 3, Nstart = 3
    real(sp), dimension(ndata)   :: xi, V, sig, Vfitted

    real(sp)                :: a(Np), F2, Fend, chisq, alamda = -0.1d0
    real(sp)                :: covar(Np,Np), alpha(Np,Np), dfit(ndata,3)
    logical                 :: mask(Np) = [.true., .true., .true.]
    
    type(pyplot)            :: plt
  
    xi = [0.0000000000000000E+00,0.6451612903225806E-01,0.1290322580645161E+00,0.1935483870967742E+00,0.2580645161290323E+00,0.3225806451612903E+00,0.3870967741935484E+00,0.4516129032258064E+00,0.5161290322580645E+00,0.5806451612903225E+00,0.6451612903225806E+00,0.7096774193548387E+00,0.7741935483870968E+00,0.8387096774193548E+00,0.9032258064516129E+00,0.9677419354838710E+00,0.1032258064516129E+01,0.1096774193548387E+01,0.1161290322580645E+01,0.1225806451612903E+01,0.1290322580645161E+01,0.1354838709677419E+01,0.1419354838709677E+01,0.1483870967741935E+01,0.1548387096774194E+01,0.1612903225806452E+01,0.1677419354838710E+01,0.1741935483870968E+01,0.1806451612903226E+01,0.1870967741935484E+01,0.1935483870967742E+01,0.2000000000000000E+01]
    V = [0.0000000000000000E+00,0.9652590349237695E-01,0.1775107606656188E+00,0.2450877415181514E+00,0.3013900160484003E+00,0.3481438130619539E+00,0.3906474563834102E+00,0.4295052157416210E+00,0.4649462049082451E+00,0.4962473370732780E+00,0.5248861954809531E+00,0.5518437021250400E+00,0.5769487178178211E+00,0.5994145652393338E+00,0.6198442180224809E+00,0.6388789909428093E+00,0.6565524894637339E+00,0.6724819427393016E+00,0.6867798604048715E+00,0.7001067911680723E+00,0.7124816543090711E+00,0.7238235436824318E+00,0.7339951036928340E+00,0.7434397486190539E+00,0.7521697712291484E+00,0.7601030676526083E+00,0.7668657897461032E+00,0.7730856669593426E+00,0.7787707791710561E+00,0.7839282458780860E+00,0.7880736011396884E+00,0.7918218312213990E+00]
    sig(1:) = 1.d1
    
    
    a(1) = (V(Nstart)-V(Nstart-1))/(xi(Nstart)-xi(Nstart-1))
    F2 = (V(Nstart+1)-V(Nstart))/(xi(Nstart+1)-xi(Nstart))
    Fend = (V(size(V))-V(size(V)-1))/(xi(size(xi))-xi(size(xi)-1))
    a(2) = abs(2./((F2-a(1))/(xi(Nstart)-xi(Nstart-1))))
    a(3) = a(1)/Fend
    
    


    call mrqmin(xi,V,sig,a,mask,covar,alpha,chisq,funcs,alamda)
    
    print *, a
    call funcs(xi,a,Vfitted,dfit)
    
    print *, a
    

    call plt%initialize(grid=.true.,xlabel='$x [nm]$',ylabel='$V (V)$', &
                            figsize=[20,15], font_size=20, title='Fitting', &
                            legend=.true.,axis_equal=.false., legend_fontsize=20, &
                            xtick_labelsize=20,ytick_labelsize=20,axes_labelsize=20)
            
    call plt%add_plot(xi,V,label='$data$', &
                    linestyle='bo',linewidth=2)
                    
    call plt%add_plot(xi,Vfitted,label='$fitting$', &
                                linestyle='r-',linewidth=2)
    call plt%savefig('fitting.png', pyfile='fitting.py')

    contains

    subroutine funcs(x,p,yfit,dydp)
        real(sp), intent(in)    :: x(:), p(:)
        real(sp), intent(out)   :: yfit(:)
        real(sp), intent(out)   :: dydp(:,:)
        
        yfit = (p(1)*p(2)*x*(p(3)-1.d0)+p(1)*x**2) / (p(3)*x+p(2)*(p(3)-1.d0))
        dydp(:,1) = (p(2) * (p(3)-1.d0)*x + x**2) / &
                    (p(2) * (p(3)-1.d0) + p(3)*x)
        dydp(:,2) = (p(1) * (x**2) * (p(3)-1.d0)**2) / &
                    (p(2) * (p(3)-1.d0) + p(3) * x)**2
        dydp(:,3) = -(p(1) * x**3) / (p(2) * (p(3)-1.d0) + p(3) * x)**2

    end subroutine funcs
    
end program testfit