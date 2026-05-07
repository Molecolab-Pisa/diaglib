module dgl_utils_c
    use dgl_interface
    use iso_c_binding

    integer, parameter :: dp = dgl_real
    integer, parameter :: ip = dgl_int

    interface check_pointer
        module procedure check_funptr
        module procedure check_ptr
    end interface

contains

    FUNCTION c_ptr_to_f_string(c_ptr_str) RESULT(f_str)
        TYPE(C_PTR), INTENT(IN) :: c_ptr_str
        CHARACTER(LEN=:), ALLOCATABLE :: f_str
        CHARACTER(KIND=C_CHAR), POINTER :: char_array(:)
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
        ALLOCATE (CHARACTER(LEN=length) :: f_str)
        DO i = 1, length
            f_str(i:i) = char_array(i)
        END DO

    END FUNCTION c_ptr_to_f_string

    subroutine dgl_error(string)
!! DiagLib error termination
        implicit none
        character(len=*), intent(in) :: string

        write (*, "(t3,a)") "-- DiagLib Error: "//string
        stop "** DiagLib issued stop signal **"
    end subroutine
!
    subroutine dgl_warning(string)
!! DiagLib error termination
        implicit none
        character(len=*), intent(in) :: string

        write (*, "(t3,a)") "-- DiagLib Warning: "//string
    end subroutine

    subroutine check_funptr(p, str)
        implicit none
        type(C_FUNPTR), intent(in) :: p
        character(len=*), intent(in) :: str

        if (.not. c_associated(p)) call dgl_error("C Pointer to "//str//" is not associated")

    end subroutine check_funptr

    subroutine check_ptr(p, str)
        implicit none
        type(C_PTR), intent(in) :: p
        character(len=*), intent(in) :: str

        if (.not. c_associated(p)) call dgl_error("C Pointer to "//str//" is not associated")

    end subroutine check_ptr

end module dgl_utils_c
