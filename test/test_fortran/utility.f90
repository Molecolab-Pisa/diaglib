module utility
implicit none

integer, parameter :: dp = selected_real_kind(15)
real(dp), parameter :: zero = 0._dp, one = 1._dp

    interface prtmat
        module procedure prtmat_r
        module procedure prtmat_i
    end interface prtmat

contains

    subroutine prtmat_r(n, m, mat)
        implicit none
        integer, intent(in) :: n, m
        real(dp), dimension(n, m), intent(in) :: mat

        character(len=5) :: s_n, s_m
        character(len=20) :: fmt
        integer :: i, j
        write (s_n, "(i0)") n
        write (s_m, "(i0)") m
        fmt = "("//s_m//"d12.3)"
        do i = 1, n
            write (*, trim(fmt)) (mat(i, j), j=1, m)
        end do
    end subroutine prtmat_r
!
    subroutine prtmat_i(n, m, mat)
        implicit none
        integer, intent(in) :: n, m
        integer, dimension(n, m), intent(in) :: mat

        character(len=5) :: s_n, s_m
        character(len=20) :: fmt
        integer :: i, j
        write (s_n, "(i0)") n
        write (s_m, "(i0)") m
        fmt = "("//s_m//"i4)"
        do i = 1, n
            write (*, trim(fmt)) (mat(i, j), j=1, m)
        end do
    end subroutine prtmat_i

end module utility