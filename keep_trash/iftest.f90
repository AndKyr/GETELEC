program iftest
    integer, parameter :: dp = 4

	real :: x(10) = [1.,2.,3.48,4.,5.32,6.,7.1,8.9,9.2,10.]
    x= [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]
    x0 = 7.36548
    
    print *, binsearch(x,x0)
    
    contains
    
    pure function binsearch(x,x0)  result(ind)
!binary search in sorted vector x for x(i) closest to x
    real(dp), intent(in)    :: x(:), x0
    integer                 :: ind, i, ia, ib, imid
       
    ia = 1
    ib = size(x)
    do i=1,size(x)
        imid = (ia + ib) / 2
        if (x(imid) < x0) then
            ia = imid
        else
            ib = imid
        endif
        if (abs(ia-ib) <= 1) then
            if (abs(x(ia)-x0) < abs(x(ib)-x0)) then
                ind = ia
            else
                ind = ib
            endif
            return
        endif
    enddo
end function binsearch

end program

