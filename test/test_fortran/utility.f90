module utility

    integer :: i, j
    integer, parameter :: dp = selected_real_kind(15)
    real(dp), parameter :: zero = 0._dp, one = 1._dp

    interface prtmat
        module procedure prtmat_r
        module procedure prtmat_i
    end interface prtmat

contains

!
    subroutine init_eigenpairs(n, n_max, eig, evec, ok)
!!
!! Make a simple guess
!!
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: n_max
        real(dp) :: eig(n_max)
        real(dp) :: evec(n, n_max)
        logical :: ok
!
        eig = zero
        evec = zero
        do i = 1, n_max
            evec(i, i) = one
        end do
        ok = .false.

    end subroutine init_eigenpairs

    subroutine dump_eigpairs(unit, n, n_targ, eig, evec, string)
!!
!! dump eigvecs to outputfile to be compares with a reference.
!!
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: n_targ
        integer, intent(in) :: unit
        character(len=*), intent(in) :: string
        real(dp), intent(inout) :: eig(n_targ), evec(n, n_targ)

        write (unit, 1000) string
        write (unit, *)
        write (unit, 1010)
        do i = 1, n_targ
            write (unit, 1020) i, eig(i)
        end do
!
! Fix the phase
!
        do i = 1, n_targ
            if (evec(1, i) .lt. zero) evec(:, i) = -evec(:, i)
        end do

        write (unit, 1031) (i, i=1, n_targ)
        do j = 1, n
            write (unit, 1022) j, (evec(j, i), i=1, n_targ)
        end do
        write (unit, *)
!
1000    format(t3, a, 1x, 'results')
1010    format(t3, 'Eigenvalues:')
1020    format(t3, i5, f14.6)
1021    format(t3, i5, f12.4)
1022    format(t3, i5, *(f14.6))
1030    format(t3, 'Eigenvector ', i3, ':')
1031    format(t3, 'Eigenvector ', tl1, *(i8, 6x))
1032    format(t3, a6, 'Eigenvector ', tl7, *(i8, 6x))
!
    end subroutine dump_eigpairs
!
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
!
end module utility
