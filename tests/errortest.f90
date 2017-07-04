program errortest

!main program designed mainly for debugging purposes
!it reads getelec emisison data from Errorfile and runs cur_dens on those data, 
! with the purpose to debug possible errors that appear with this specific data set.

use GeTElEC, only: EmissionData, dp, cur_dens, print_data, plot_barrier, debug

implicit none

type(EmissionData)      :: this
real(dp)                :: x(512), V(512)
character               :: dash
character(len=32)       :: line, str, reg_str, sharp_str, approx_str
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
read (fid,'(/A15,A20,/A15,A20,/A15,I20,/A15,A20,/A15,I20)') str, reg_str, str, &
        sharp_str,  str, this%mode, str, approx_str,    str, this%ierr

read (fid,'(A32)') line

if ( line(1:2) /= '--') then
    read (fid,'(A32)') line
    do i = 1, 512
        read (fid,'(A32)') line
        if ( line(1:2) == '--') exit
        read(line,'(i3,ES14.4,ES15.5)') j, x(i), V(i)
    enddo
endif

close(fid)


if (trim(reg_str) == 'Field') this%regime = 1 
if (trim(reg_str) == 'Intermediate') this%regime = 0
if (trim(reg_str) == 'Thermal') this%regime = -1

if (trim(sharp_str) == 'Sharp') this%sharpness = 1 
if (trim(sharp_str) == 'Blunt') this%sharpness = 0

if (trim(approx_str) == 'Richardson-Dushman') this%approx = -2
if (trim(approx_str) == 'Miller-Good') this%approx = -1
if (trim(approx_str) == 'General T-F (Jensen)') this%approx = 0
if (trim(approx_str) == 'Automatic selection') this%approx = 1
if (trim(approx_str) == 'Full integration') this%approx = 2

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
