program errortest

!main program designed mainly for debugging purposes
!it reads getelec emisison data from Errorfile and runs cur_dens on those data, 
! with the purpose to debug possible errors that appear with this specific data set.

use GeTElEC, only: EmissionData, dp, cur_dens, print_data, plot_barrier, debug

implicit none

type(EmissionData)      :: this
real(dp)                :: x(512), V(512)
character               :: dash
character(len=32)       :: line, str
integer                 :: fid, i,j, Nr

    open(fid, file = 'errors.txt', action = 'read')
    
    read (fid,'(A32)') line
    read (fid,'(A10,ES12.4,A10/,A10,ES12.4,A10,/A10,ES12.4)') str, this%F, str, &
            str, this%R, str, str, this%gamma
    read (fid,'(A10,ES12.4,A10,/A10,ES12.4,A10)') str, this%W, str,str, this%kT
    read (fid,'(/A10,ES12.4,A10,/A10,ES12.4,A10,/A10,A10,A10)') str, this%Jem, str, &
            str, this%heat, str, str, str, str
    read (fid,'(/A10,ES12.4,A10,/A10,ES12.4,A10,/A10,ES12.4,A10)') str, this%xm, &
            str, str, this%xmax, str, str,  this%Um, str
    read (fid,'(A10,ES12.4,/A10,ES12.4,A10,/A10,ES12.4,A10)') str, this%Gam, str, &
            this%maxbeta, str, str, this%minbeta, str               
    read (fid,'(//A21,A1,/A21,A1,/A10,I12,/A10,I12,/A10,I12)') str, this%regime, &
            str, this%sharpness, str, this%mode, str, this%approx, str, this%ierr
      
    read(fid, '(A32)')  str
    
    do i = 1, 512
        read (fid,'(A32)') line
        if ( line(1:2) == '--') exit
        read(line,'(i3,ES14.4,ES15.5)') j, x(i), V(i)
    enddo
    
    close(fid)
    
    Nr = i-1
    allocate(this%xr(Nr),this%Vr(Nr))
    this%xr = x(:Nr)
    this%Vr = V(:Nr)
    
    print *, 'Read the following emission data:'
    print *, 'Nr = ', Nr
    call print_data(this,.true.)
    
    call cur_dens(this)
    
    print *, 'After recalculation:'
    
    call print_data(this)


end program errortest
