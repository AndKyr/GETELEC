program typetest

type    :: malakas
    real, allocatable      :: x(:),y(:)
end type malakas

type(malakas)   :: mal

allocate(mal%x(3), mal%y(3))

mal%x = [1.,2.,3.]
mal%y = [1.,2.,3.]

print *, mal%x

end program typetest
