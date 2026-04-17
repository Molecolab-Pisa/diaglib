module dgl_utils_c
    use dgl_interface
    use iso_c_binding

    integer, parameter :: dp = dgl_real

contains

FUNCTION c_ptr_to_f_string(c_ptr_str) RESULT(f_str)
        TYPE(C_PTR), INTENT(IN)           :: c_ptr_str
        CHARACTER(LEN=:), ALLOCATABLE     :: f_str
        CHARACTER(KIND=C_CHAR), POINTER   :: char_array(:)
        INTEGER :: length, i

        ! Step 1: Check for NULL pointer
        IF (.NOT. C_ASSOCIATED(c_ptr_str)) THEN
            f_str = ''
            RETURN
        END IF

        ! Step 2: Associate C pointer with Fortran pointer
        CALL C_F_POINTER(c_ptr_str, char_array, [HUGE(0)])

        ! Step 3: Scan for null terminator to get length
        length = 0
        DO WHILE (char_array(length + 1) /= C_NULL_CHAR)
            length = length + 1
        END DO

        ! Step 4: Copy characters into allocatable Fortran string
        ALLOCATE(CHARACTER(LEN=length) :: f_str)
        DO i = 1, length
            f_str(i:i) = char_array(i)
        END DO

END FUNCTION c_ptr_to_f_string

end module dgl_utils_c