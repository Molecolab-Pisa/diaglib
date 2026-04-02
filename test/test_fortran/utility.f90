module utility
!!
!! Common module for testing
!!
    integer, parameter :: dp = selected_real_kind(15)
    real(dp), parameter :: zero = 0._dp, one = 1._dp, two = 2._dp, five = 5._dp

    integer, parameter :: n = 500, n_targ = 5, n_max = 10
    integer, parameter :: max_iter = 100, dav_iter = 10
    logical :: verbose = .false.
    real(dp), parameter :: tol = 1.0e-10_dp, shift = 0.0_dp
    integer, parameter :: memory = 100
    character(len=2), parameter :: memory_unit = "MB"
    procedure(), pointer :: mx_p => null()

    character(len=30) :: reference_fname = "output_reference.txt"
    character(len=30) :: results_fname = "output_fortran.txt"
    integer, parameter :: lutest = 100, luref = 101
    character(len=20), parameter :: f_string = "('-- ', a, /)"
    integer :: i, j

    interface prtmat
        module procedure prtmat_r
        module procedure prtmat_i
    end interface prtmat

contains

!
    subroutine init_eigenpairs(ld, n_eig, eig, evec, ok)
!!
!! Make a simple guess
!!
        implicit none
        integer, intent(in) :: ld
        integer, intent(in) :: n_eig
        real(dp) :: eig(n_eig)
        real(dp) :: evec(ld, n_eig)
        logical :: ok
!
        eig = zero
        evec = zero
        do i = 1, n_eig
            evec(i, i) = one
        end do
        ok = .false.

    end subroutine init_eigenpairs

    logical function compare_eigs(ld, n_eig, thresh, eig, evec, string) result(success)
        integer, intent(in) :: ld, n_eig
        real(dp), intent(in) :: thresh
        real(dp), intent(in) :: eig(n_eig)
        real(dp), intent(inout) :: evec(ld, n_eig)
        character(len=*) :: string

        real(dp), allocatable :: ex_eig(:), ex_evec(:, :), norms(:, :)

        success = .true.

        allocate (ex_eig(n_eig), ex_evec(ld, n_eig), norms(n_eig, 2))

        call read_reference(ld, n_eig, ex_eig, ex_evec, string)

        if (abs(maxval(eig - ex_eig)) .gt. thresh) success = .false.

        do i = 1, n_eig
            if (evec(1, i) .lt. zero) evec(:, i) = -evec(:, i)
            norms(i, 1) = sqrt(dot_product(evec(:, i), evec(:, i)))
            norms(i, 2) = sqrt(dot_product(evec(:, i), ex_evec(:, i)))
        end do
        if (abs(maxval(norms(:, 1) - norms(:, 2))) .gt. thresh) success = .false.

        if (.not. success) then
            write (*, "(t3,a,*(d14.4))") "Computed  Eigenvals:", eig
            write (*, "(t3,a,*(d14.4))") "Reference Eigenvals:", ex_eig
            write (*, "(t3,a,*(d14.4))") "Difference:         ", abs(ex_eig - eig)
            write (*, *)
            write (*, "(t3,a,*(d14.4))") "Norm of Computed  Eigenvecs:     ", norms(:, 1)
            write (*, "(t3,a,*(d14.4))") "Overlap with Reference Eigenvecs:", norms(:, 2)
            write (*, "(t3,a,*(d14.4))") "Difference:                      ", abs(norms(:, 2) - norms(:, 1))
        end if

        deallocate (ex_eig, ex_evec, norms)

    end function compare_eigs

    subroutine check_reference(ld, n_eig)
        implicit none
        integer, intent(in) :: ld, n_eig

        integer :: n_eig_read, ld_read
        character(len=200) :: line

        if (.not. file_exist(trim(reference_fname))) then
            write (*, "(t3,a)") "Reference file does not exist, first run the 'reference' executable"
            stop
        end if

        open (file=trim(reference_fname), unit=luref, status="old")

        do
            read (luref, "(a)") line
            if (index(line, "Eigenvalues") .ne. 0) exit
        end do

        n_eig_read = 0
        do
            read (luref, "(a)") line
            if (index(line, "Eigenvectors") .ne. 0) exit
            n_eig_read = n_eig_read + 1
        end do

        if (n_eig .ne. n_eig_read) then
            write (*, "(t3,a)") "Reference file contains the wrong number of eigenvalues, "// &
                "re-run the 'reference' executable"
            stop
        end if

        ld_read = 0
        do
            read (luref, "(a)") line
            if (len_trim(line) .eq. 0) exit
            ld_read = ld_read + 1
        end do

        if (ld .ne. ld_read) then
            write (*, "(t3,a)") "Reference file contains results for a differently sized matrix, "// &
                "re-run the 'reference' executable"
            stop
        end if

        close (luref)

    end subroutine check_reference

    subroutine read_reference(ld, n_eig, eig, evec, string)
        integer, intent(in) :: ld, n_eig
        real(dp), intent(inout) :: eig(n_eig), evec(ld, n_eig)
        character(len=*), intent(in) :: string

        character(len=200) :: line
        integer :: k

        open (file=trim(reference_fname), unit=luref, status="old")

        do
            read (luref, "(a)") line
            if (index(line, string) .ne. 0) exit
        end do

        read (luref, *)
        read (luref, *)
        do i = 1, n_eig
            read (luref, *) k, eig(i)
        end do

        read (luref, *)
        do i = 1, ld
            read (luref, *) k, (evec(i, j), j=1, n_eig)
        end do

        close (luref)

    end subroutine read_reference

    subroutine sort_eigenpairs(ld, n_eig, eig, evec)
        implicit none
        integer, intent(in) :: ld
        integer, intent(in) :: n_eig
        real(dp), intent(inout) :: eig(ld)
        real(dp), intent(inout) :: evec(ld, ld, 2)

        logical, allocatable :: mask(:)
        integer :: idx(1)
        real(dp) :: copy
        real(dp), allocatable :: copy_arr(:)

        allocate (mask(ld), copy_arr(ld))
        mask = .true.

        do i = 1, n_eig
            idx = minloc(eig, mask=mask)
            copy = eig(i)
            eig(i) = eig(idx(1))
            eig(idx(1)) = copy
            mask(i) = .false.

            copy_arr = evec(:, i, 1)
            evec(:, i, 1) = evec(:, idx(1), 1)
            evec(:, idx(1), 1) = copy_arr

            copy_arr = evec(:, i, 2)
            evec(:, i, 2) = evec(:, idx(1), 2)
            evec(:, idx(1), 2) = copy_arr
        end do

        deallocate (mask, copy_arr)

    end subroutine sort_eigenpairs

    subroutine dump_eigpairs(unit, ld, n_eig, eig, evec, string)
!!
!! dump eigvecs to outputfile to be compares with a reference.
!!
        implicit none
        integer, intent(in) :: ld
        integer, intent(in) :: n_eig
        integer, intent(in) :: unit
        character(len=*), intent(in) :: string
        real(dp), intent(inout) :: eig(n_eig), evec(ld, n_eig)

        if (.not. file_opened(unit)) then
            write (*, "(t3,a)") "Tried to dump on file that is not open"
            stop
        end if

        write (unit, 1000) string
        write (unit, *)
        write (unit, 1010)
        do i = 1, n_eig
            write (unit, 1020) i, eig(i)
        end do
!
! Fix the phase
!
        do i = 1, n_eig
            if (evec(1, i) .lt. zero) evec(:, i) = -evec(:, i)
        end do

        write (unit, 1030) (i, i=1, n_eig)
        do j = 1, ld
            write (unit, 1022) j, (evec(j, i), i=1, n_eig)
        end do
        write (unit, *)
!
1000    format(/, t3, a, ' results')
1010    format(t3, 'Eigenvalues')
1020    format(t3, i7, f22.16)
1022    format(t3, i7, *(f22.16))
1030    format(t3, 'Eigenvectors ', *(i16, 6x))
!
    end subroutine dump_eigpairs

    logical function file_exist(fname)
        implicit none
        character(len=*), intent(in) :: fname

        inquire (file=fname, exist=file_exist)

    end function file_exist
!
    logical function file_opened(unit)
        implicit none
        integer, intent(in) :: unit

        inquire (unit=unit, opened=file_opened)

    end function file_opened

    subroutine prtmat_r(ld, m, mat)
        implicit none
        integer, intent(in) :: ld, m
        real(dp), dimension(ld, m), intent(in) :: mat

        character(len=5) :: s_n, s_m
        character(len=20) :: fmt

        write (s_n, "(i0)") ld
        write (s_m, "(i0)") m
        fmt = "("//s_m//"d12.3)"
        do i = 1, ld
            write (*, trim(fmt)) (mat(i, j), j=1, m)
        end do

    end subroutine prtmat_r
!
    subroutine prtmat_i(ld, m, mat)
        implicit none
        integer, intent(in) :: ld, m
        integer, dimension(ld, m), intent(in) :: mat

        character(len=5) :: s_n, s_m
        character(len=20) :: fmt

        write (s_n, "(i0)") ld
        write (s_m, "(i0)") m
        fmt = "("//s_m//"i4)"
        do i = 1, ld
            write (*, trim(fmt)) (mat(i, j), j=1, m)
        end do

    end subroutine prtmat_i
!
end module utility
