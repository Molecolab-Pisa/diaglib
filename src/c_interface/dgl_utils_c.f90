module dgl_utils_c
    use dgl_interface
    use iso_c_binding
    implicit none

    integer, parameter :: dp = dgl_real
    integer, parameter :: ip = dgl_int
!
! kind of the integers exchanged with C: a fixed-width type that matches ip
! and the dgl_int typedef in diaglib.h
!
#ifdef DGL_INT_KIND_8
    integer, parameter :: c_ip = C_INT64_T
#elif DGL_INT_KIND_4
    integer, parameter :: c_ip = C_INT32_T
#endif

    interface pointer_ok
        module procedure funptr_ok
        module procedure ptr_ok
    end interface

!
! interfaces of the C routines supplied by the user, as declared in diaglib.h
!
    abstract interface
        subroutine c_matvec(n, m, x, y) bind(C)
            import :: c_ip, C_DOUBLE
            integer(c_ip), intent(in) :: n, m
            real(C_DOUBLE), intent(in) :: x(n, m)
            real(C_DOUBLE), intent(inout) :: y(n, m)
        end subroutine c_matvec
!
        subroutine c_precnd(n, m, shift, x, y) bind(C)
            import :: c_ip, C_DOUBLE
            integer(c_ip), intent(in) :: n, m
            real(C_DOUBLE), intent(in) :: shift
            real(C_DOUBLE), intent(in) :: x(n, m)
            real(C_DOUBLE), intent(inout) :: y(n, m)
        end subroutine c_precnd
!
        subroutine c_lrprec(n, m, fac, xp, xm, yp, ym) bind(C)
            import :: c_ip, C_DOUBLE
            integer(c_ip), intent(in) :: n, m
            real(C_DOUBLE), intent(in) :: fac
            real(C_DOUBLE), intent(in) :: xp(n, m), xm(n, m)
            real(C_DOUBLE), intent(inout) :: yp(n, m), ym(n, m)
        end subroutine c_lrprec
    end interface

contains

    function integer_kind_c() result(kind_bytes) bind(C, name="dgl_integer_kind")
!! Size in bytes of the integers of this DiagLib build (4 or 8), so that the users of the
!! compiled library (e.g., the python interface) do not have to know how it was built.
        implicit none
        integer(C_INT) :: kind_bytes
        kind_bytes = int(storage_size(0_c_ip)/8, C_INT)
    end function integer_kind_c

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

    logical function funptr_ok(p, str)
!! Check that a C function pointer is associated, printing an error if it is not
        implicit none
        type(C_FUNPTR), intent(in) :: p
        character(len=*), intent(in) :: str

        funptr_ok = c_associated(p)
        if (.not. funptr_ok) write (*, "(t3,a)") "-- DiagLib Error: C Pointer to "//str//" is not associated"

    end function funptr_ok

    logical function ptr_ok(p, str)
!! Check that a C pointer is associated, printing an error if it is not
        implicit none
        type(C_PTR), intent(in) :: p
        character(len=*), intent(in) :: str

        ptr_ok = c_associated(p)
        if (.not. ptr_ok) write (*, "(t3,a)") "-- DiagLib Error: C Pointer to "//str//" is not associated"

    end function ptr_ok

end module dgl_utils_c
