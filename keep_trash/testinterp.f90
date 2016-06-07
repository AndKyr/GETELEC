program test

integer, parameter  :: dp = 8    

real(dp)                            :: x , xi(32), yi(32), y(2)

xi = linspace(0.d0,1.d0,32)
yi =  xi**2 + 3.d0 * xi
x = 0.1598752
y = interp1(xi,yi,x)
print *, y
print *,   x**2 + 3.d0 * x

contains

pure function interp1(xi,yi,x) result(y)
!simple linear interpolation function (same as interp1 of matlab)
!appropriate for non-uniform linspace

    real(dp), intent(in)    :: xi(:), yi(:), x
    real(dp)                :: y(2) ! 1: output, 2: flag
    integer                 :: dN, Nclose, binout(2)
    
    binout = binsearch(xi,x)
    y(2) = real(binout(2),dp)
    if (binout(2) /= 0) return
    dN = 1
    Nclose = binout(1)
    if (xi(Nclose) > x) dN = -1
    y(1) = yi(Nclose) + (x-xi(Nclose)) * &
        ((yi(Nclose + dN) - yi(Nclose)) / (xi(Nclose + dN) - xi(Nclose)))
    
end function interp1

pure function binsearch(x, x0)  result(ind)
!binary search in sorted vector x for x(i) closest to x
    real(dp), intent(in)    :: x(:), x0
    integer                 :: ind(2) ! 1: output, 2: iflag
    !iflag 0: fine, 1: x0 out of bounds, -1: x not sorted
    integer                 :: i, ia, ib, imid
    
    if ((x0 < x(1)) .or. (x0 > x(size(x)))) then
        ind(2) = 1
        return
    endif
    
    ia = 1
    ib = size(x) +1
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
            ind(2) = 0
            return
        endif
    enddo
    ind(2) = -1
    
end function binsearch

pure function linspace(a,b,N) result(x)
    real(dp), intent(in)    ::a,b
    integer, intent(in)     ::N
    real(dp)                :: dx,x(N)
    integer                 :: i
    
    dx=(b-a)/(N-1)
    do i=1,N
        x(i)=a+(i-1)*dx
    enddo
end function linspace
end program










	
