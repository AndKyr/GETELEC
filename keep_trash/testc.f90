module skata

use iso_c_binding

type, bind(c)   :: pass
    integer(c_int)  :: lenc, lenf
    type(c_ptr)     :: c, f
end type pass

contains

subroutine foo(arrays) bind(c)

    
    type(pass), intent(in)  :: arrays
    real(c_double), pointer :: c_array(:)
    
    call c_f_pointer(arrays%c, c_array, [arrays%lenc])
    
    print *, c_array
end subroutine foo

end module skata
