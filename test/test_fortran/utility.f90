module utility
!!
!! Common module for testing
!!
#ifdef DGL_INT_KIND_8
    integer, parameter :: ip = selected_int_kind(15)
#elif DGL_INT_KIND_4
    integer, parameter :: ip = selected_int_kind(8)
#endif
    integer(ip), parameter :: dp = selected_real_kind(15)
    real(dp), parameter :: zero = 0._dp, one = 1._dp, two = 2._dp, five = 5._dp, half = 0.5_dp

    integer(ip), parameter :: n = 700, n_targ = 5, n_max = 10
    integer(ip), parameter :: max_iter = 100, dav_iter = 10
    logical :: verbose = .false.
    real(dp), parameter :: tol = 1.0e-10_dp, shift = 0.0_dp
!
! thresholds used to compare the results with the reference: maximum error on the eigenvalues
! and maximum distance between the normalized computed and reference eigenvectors (about the
! sine of the angle between them). with tol = 1e-10, the observed errors are at most about
! 2e-12 on the eigenvalues and 6e-11 on the eigenvectors, for all drivers and guesses.
!
    real(dp), parameter :: eig_thresh = 1.0e-9_dp, vec_thresh = 1.0e-7_dp
!
! number of runs for each test that starts from a random guess
!
    integer(ip), parameter :: n_random_runs = 3
    integer(ip), parameter :: memory = 100
    character(len=2), parameter :: memory_unit = "MB"
    procedure(), pointer :: mx_p => null()

    character(len=30) :: reference_fname = "output_reference.txt"
    character(len=30) :: results_fname = "output_fortran.txt"
    integer(ip), parameter :: lutest = 100, luref = 101
    character(len=20), parameter :: f_string = "('-- ', a, /)"
    integer(ip) :: i, j

    interface prtmat
        module procedure prtmat_r
        module procedure prtmat_i
    end interface prtmat

contains

!
    subroutine init_eigenpairs(ld, n_eig, eig, evec, ok, random_guess)
!!
!! Make a simple guess (unit vectors), or no guess at all (zero vectors), in which case
!! DiagLib starts from a random guess
!!
        implicit none
        integer(ip), intent(in) :: ld
        integer(ip), intent(in) :: n_eig
        real(dp) :: eig(n_eig)
        real(dp) :: evec(ld, n_eig)
        logical :: ok
        logical, intent(in) :: random_guess
!
        integer(ip) :: k
!
        eig = zero
        evec = zero
        if (.not. random_guess) then
            do k = 1, n_eig
                evec(k, k) = one
            end do
        end if
        ok = .false.

    end subroutine init_eigenpairs

    real(dp) function close_matrix_element(ii, jj) result(val)
!!
!! Element (ii,jj) of a non-symmetric matrix with a real spectrum and a nearly degenerate
!! pair of low eigenvalues (about 2.950 and 3.056), used to test the root following of the
!! non-symmetric Davidson driver. Used both to build the reference and in the matvecs.
!!
        implicit none
        integer(ip), intent(in) :: ii, jj
        real(dp), parameter :: diag(5) = [2.0_dp, 3.0_dp, 3.001_dp, 3.5_dp, 4.2_dp]
!
        if (ii .eq. jj) then
            if (ii .le. 5) then
                val = diag(ii)
            else
                val = real(ii + 1, dp)
            end if
        else
            val = 0.3_dp*(real(ii, dp)/real(jj, dp))/real(ii + jj, dp)
            if (mod(ii + jj - 2_ip, 3_ip) .eq. 0) val = -val
        end if

    end function close_matrix_element

    real(dp) function cplx_matrix_element(ii, jj) result(val)
!!
!! Element (ii,jj) of the non-symmetric test matrix applied by arx, plus an antisymmetric
!! coupling between the second and third basis functions, which produces a pair of complex
!! eigenvalues (about 3.49 +- 0.84 i) among the lowest real ones. The drivers have to skip them.
!!
        implicit none
        integer(ip), intent(in) :: ii, jj
!
        if (ii .eq. jj) then
            val = real(ii + 1, dp)
        else if (ii .eq. 2 .and. jj .eq. 3) then
            val = one
        else if (ii .eq. 3 .and. jj .eq. 2) then
            val = -one
        else
            val = (real(ii, dp)/real(jj, dp))/real(ii + jj, dp)
        end if

    end function cplx_matrix_element

    logical function compare_eigs(ld, n_eig, eig, evec, string) result(success)
!!
!! Compare computed eigenpairs with the reference ones. The eigenvectors are compared
!! through the distance between the normalized vectors (with the sign that minimizes it),
!! which does not depend on their normalization and phase.
!!
        integer(ip), intent(in) :: ld, n_eig
        real(dp), intent(in) :: eig(n_eig)
        real(dp), intent(in) :: evec(ld, n_eig)
        character(len=*) :: string

        real(dp), allocatable :: ex_eig(:), ex_evec(:, :), u(:), v(:)
        real(dp) :: eig_err, vec_err, dist
        integer(ip) :: k

        allocate (ex_eig(n_eig), ex_evec(ld, n_eig), u(ld), v(ld))

        call read_reference(ld, n_eig, ex_eig, ex_evec, string)

        eig_err = maxval(abs(eig - ex_eig))

        vec_err = zero
        do k = 1, n_eig
            u = evec(:, k)/sqrt(dot_product(evec(:, k), evec(:, k)))
            v = ex_evec(:, k)/sqrt(dot_product(ex_evec(:, k), ex_evec(:, k)))
            dist = min(sqrt(dot_product(u - v, u - v)), sqrt(dot_product(u + v, u + v)))
            vec_err = max(vec_err, dist)
        end do
!
! written this way, a NaN makes the comparison fail
!
        success = eig_err .le. eig_thresh .and. vec_err .le. vec_thresh

        write (*, "(t3,a,d10.2,a,d10.2)") "max eigenvalue error:", eig_err, "   max eigenvector error:", vec_err
        if (.not. success) then
            write (*, "(t3,a,*(d14.4))") "Computed  Eigenvals:", eig
            write (*, "(t3,a,*(d14.4))") "Reference Eigenvals:", ex_eig
            write (*, "(t3,a,*(d14.4))") "Difference:         ", abs(ex_eig - eig)
        end if

        deallocate (ex_eig, ex_evec, u, v)

    end function compare_eigs

    subroutine check_reference(ld)
        implicit none
        integer(ip), intent(in) :: ld

        integer(ip) :: ld_read
        character(len=200) :: line

        if (.not. file_exist(trim(reference_fname))) then
            write (*, "(t3,a)") "Reference file does not exist, first run the 'reference' executable"
            stop 1
        end if

        open (file=trim(reference_fname), unit=luref, status="old")

        do
            read (luref, "(a)") line
            if (index(line, "Eigenvalues") .ne. 0) exit
        end do

        ld_read = 0
        do
            read (luref, "(a)") line
            if (index(line, "Eigenvectors") .ne. 0) exit
            ld_read = ld_read + 1
        end do

        if (ld .ne. ld_read) then
            write (*, "(t3,a)") "Reference file contains results for a differently sized matrix, "// &
                "re-run the 'reference' executable"
            stop 1
        end if

        close (luref)

    end subroutine check_reference

    subroutine read_reference(ld, n_eig, eig, evec, string)
        integer(ip), intent(in) :: ld, n_eig
        real(dp), intent(inout) :: eig(n_eig), evec(ld, n_eig)
        character(len=*), intent(in) :: string

        character(len=200) :: line
        integer(ip) :: k

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

        do
            read (luref, "(a)") line
            if (index(line, "Eigenvectors") .ne. 0) exit
        end do

        do i = 1, ld
            read (luref, *) k, (evec(i, j), j=1, n_eig)
        end do

        close (luref)

    end subroutine read_reference

    subroutine sort_eigenpairs(ld, n_eig, eig, evec)
!!
!! sort the first n_eig eigenvalues in ascending order, and the two sets of
!! eigenvectors with them.
!!
        implicit none
        integer(ip), intent(in) :: ld
        integer(ip), intent(in) :: n_eig
        real(dp), intent(inout) :: eig(ld)
        real(dp), intent(inout) :: evec(ld, ld, 2)

        integer(ip) :: k, l, i_min
        real(dp) :: copy
        real(dp), allocatable :: copy_arr(:)

        allocate (copy_arr(ld))
!
! plain selection sort, with local loop indices. The same algorithm used to be
! written with minloc(eig, mask=mask, kind=ip), which is not usable: with both a
! mask and a kind argument, gfortran calls a libgfortran routine that mishandles
! its back argument and returns the LAST of several equal minima, where the
! standard requires the first (checked with gfortran 11.4 and 12.3 at -O0; from
! -O1 the call is inlined and correct, which is why only Debug builds were
! affected). On some builds of libgfortran the same call aborts instead, with
! "Assertion 'back == 0' failed", which took down the whole reference generator.
!
        do k = 1, n_eig
            i_min = k
            do l = k + 1, ld
                if (eig(l) .lt. eig(i_min)) i_min = l
            end do
            if (i_min .eq. k) cycle

            copy = eig(k)
            eig(k) = eig(i_min)
            eig(i_min) = copy

            copy_arr = evec(:, k, 1)
            evec(:, k, 1) = evec(:, i_min, 1)
            evec(:, i_min, 1) = copy_arr

            copy_arr = evec(:, k, 2)
            evec(:, k, 2) = evec(:, i_min, 2)
            evec(:, i_min, 2) = copy_arr
        end do

        deallocate (copy_arr)

    end subroutine sort_eigenpairs

    subroutine dump_eigpairs(unit, ld, n_eig, eig, evec, string)
!!
!! dump eigvecs to outputfile to be compares with a reference.
!!
        implicit none
        integer(ip), intent(in) :: ld
        integer(ip), intent(in) :: n_eig
        integer(ip), intent(in) :: unit
        character(len=*), intent(in) :: string
        real(dp), intent(inout) :: eig(n_eig), evec(ld, n_eig)

        if (.not. file_opened(unit)) then
            write (*, "(t3,a)") "Tried to dump on file that is not open"
            stop 1
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
        integer(ip), intent(in) :: unit

        inquire (unit=unit, opened=file_opened)

    end function file_opened

    subroutine prtmat_r(ld, m, mat)
        implicit none
        integer(ip), intent(in) :: ld, m
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
        integer(ip), intent(in) :: ld, m
        integer(ip), dimension(ld, m), intent(in) :: mat

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
