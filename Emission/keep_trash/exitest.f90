program exitest

integer i

do i=1,10
    print *, i
    if (i==8) exit
enddo

print *,i

end program exitest