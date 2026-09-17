program test_fortran
!!
!! tests the various functionalities of diaglib, comparing the results with the reference
!! computed with lapack by the dgl_reference program.
!!
!! 1. all the drivers are run from a simple guess (unit vectors) and n_random_runs times from a
!!    random guess (zero vectors in input, so that diaglib generates the guess). results obtained
!!    from a random guess are not reproducible, so only the converged eigenpairs are checked.
!!    the non-symmetric driver is run for right, left and both eigenvectors, and for three
!!    matrices: one with well separated eigenvalues, one with nearly degenerate eigenvalues, in
!!    which the roots change order during the iterations and have to be followed, and one with
!!    a pair of complex eigenvalues among the lowest ones, which have to be skipped.
!! 2. all the drivers are run with a level shift, in verbose mode, from a guess that is not
!!    orthonormal (LOBPCG with a shifted preconditioner, as the Davidson drivers), and from an
!!    incomplete guess (a duplicated vector and zero vectors).
!! 3. the non-symmetric driver is run with a poor preconditioner, so that the expansion space
!!    is restarted several times, also when computing the left eigenvectors after the right ones.
!! 4. invalid input (including NaN returned by the user-supplied routines) has to be reported
!!    through dgl_info, and values of dav_iter that are too small or too large have to be handled.
!! 5. DiagLib is called from inside a matrix-vector callback: the nested calls must not
!!    interfere with the outer one.
!! 6. all the drivers that use (b_)ortho_vs_x are run on a tridiagonal matrix, starting from
!!    unit vectors, so that the preconditioned residuals are linearly dependent.
!! 7. Davidson and LOBPCG are run on a matrix with known eigenvalues whose eigenvectors are far
!!    from unit vectors (a diagonal matrix rotated by a Householder reflection), so that the
!!    diagonal preconditioner is only approximate, with and without the shift of the
!!    preconditioner.
!!
    use dgl_interface
    use direct_matvecs
    use utility
    implicit none
!
! kinds of guess
!
    integer(ip), parameter :: simple_guess = 0, random_guess = 1, scaled_guess = 2, partial_guess = 3
!
! non-symmetric test matrices
!
    integer(ip), parameter :: separated = 1, close_eigs = 2, cplx_eigs = 3
!
! options of the current tests
!
    integer(ip) :: guess_kind
    real(dp) :: test_shift
    logical :: test_verbose, weak_precnd
!
! LOBPCG is run with the default (no shift for the preconditioner), except in pass 2
!
    logical :: lobpcg_precnd_shift = .false.
!
    logical :: ok
    integer(ip) :: n_tests = 0, n_failed = 0, run, i_side, i_mat
    real(dp), allocatable :: eig(:), evec(:, :), evec_2(:, :)
    character(len=2), parameter :: sides(3) = ["R ", "L ", "LR"]
!
! Check if the reference exists
!
    call check_reference(n)
!
! open a text file to dump the output.
!
    open (file=trim(results_fname), form='formatted', access='sequential', &
          unit=lutest, status='unknown')
!
! allocate space for eigenvalues and eigenvectors
!
    allocate (eig(n_max), evec(n, n_max), evec_2(n, n_max))
!
! Associate pointer for generalized problems to metric-vector
! product routine
!
    mx_p => mx
!
!====================================
! TESTING AREA
!====================================
!
! 1. simple and random guesses
!
    test_shift = shift
    test_verbose = verbose
    weak_precnd = .false.
    do run = 0, n_random_runs
        guess_kind = merge(random_guess, simple_guess, run .gt. 0)
        call run_all_drivers()
    end do
!
! 2. level shift, verbose output, guess that is not orthonormal
!
    guess_kind = scaled_guess
    test_shift = 0.5_dp
    test_verbose = .true.
    lobpcg_precnd_shift = .true.
    call run_all_drivers()
    lobpcg_precnd_shift = .false.
    guess_kind = partial_guess
    test_shift = shift
    test_verbose = verbose
    call run_all_drivers()
!
! 3. poor preconditioner, many restarts
!
    guess_kind = simple_guess
    test_shift = shift
    test_verbose = verbose
    weak_precnd = .true.
    do i_side = 1, size(sides)
        call test_nonsym_davidson(sides(i_side), close_eigs)
    end do
    weak_precnd = .false.
!
! 4. invalid input and extreme values of dav_iter
!
    call test_input_errors()
    call test_dav_iter()
!
! 5. nested calls
!
    call test_nested()
!
! 6. linearly dependent preconditioned residuals
!
    call test_dependent_residuals()
!
! 7. rotated matrix, approximate preconditioner
!
    call test_rotated()
!
! close the output file:
!
    close (lutest)
!
! free the memory:
!
    deallocate (evec, evec_2, eig)
!
    write (6, "(t3,a,i0,a,i0,a)") "Summary: ", n_failed, " failed tests out of ", n_tests, "."
    if (n_failed .gt. 0) stop 1
!
contains
!
    subroutine run_all_drivers()
        implicit none
        call test_davidson(.false.)
        call test_davidson(.true.)
        call test_lobpcg(.false.)
        call test_lobpcg(.true.)
        do i_mat = separated, cplx_eigs
            do i_side = 1, size(sides)
                call test_nonsym_davidson(sides(i_side), i_mat)
            end do
        end do
        call test_smogd()
    end subroutine run_all_drivers
!
    subroutine init_guess(ld, vecs)
!
! build the guess for the current test
!
        implicit none
        integer(ip), intent(in) :: ld
        real(dp), intent(inout) :: vecs(ld, n_max)
        integer(ip) :: k
!
        call init_eigenpairs(ld, n_max, eig, vecs, ok, guess_kind .eq. random_guess)
        if (guess_kind .eq. scaled_guess) then
            do k = 1, n_max
                vecs(k, k) = 2.0_dp
                vecs(k + 1, k) = 0.5_dp
            end do
        end if
!
! guess only for the first n_targ vectors: then a copy of the first one, and zero vectors
!
        if (guess_kind .eq. partial_guess) then
            vecs(:, n_targ + 1) = vecs(:, 1)
            vecs(:, n_targ + 2:) = zero
        end if
    end subroutine init_guess
!
    function test_label(driver) result(label)
!
! describe the current test
!
        implicit none
        character(len=*), intent(in) :: driver
        character(len=120) :: label
!
        select case (guess_kind)
        case (random_guess)
            label = driver//", random guess"
        case (scaled_guess)
            label = driver//", shifted, verbose, non-orthonormal guess"
        case (partial_guess)
            label = driver//", incomplete guess"
        case default
            label = driver//", simple guess"
        end select
        if (weak_precnd) label = trim(label)//", poor preconditioner"
    end function test_label
!
    subroutine check_result(label, ok, info, compared)
!
! record the result of a test: the driver has to converge without errors,
! and the eigenpairs have to agree with the reference
!
        implicit none
        character(len=*), intent(in) :: label
        logical, intent(in) :: ok, compared
        integer(ip), intent(in) :: info
!
        n_tests = n_tests + 1
        if (info .eq. dgl_success .and. ok .and. compared) then
            write (6, f_string) label//": PASSED"
        else
            n_failed = n_failed + 1
            if (info .ne. dgl_success) then
                write (6, "('-- ', a, ': FAILED (error ', i0, ')', /)") label, info
            else if (.not. ok) then
                write (6, f_string) label//": FAILED (not converged)"
            else
                write (6, f_string) label//": FAILED (results do not match the reference)"
            end if
        end if
    end subroutine check_result
!
    subroutine test_davidson(generalized)
        implicit none
        logical, intent(in) :: generalized
        integer(ip) :: info
        logical :: same
        character(len=120) :: label
!
        if (generalized) then
            label = test_label("Generalized Davidson")
        else
            label = test_label("Davidson")
        end if
        call init_guess(n, evec)
!
        write (6, f_string) 'testing '//trim(label)//':'
        if (generalized) then
            call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
                                     dgl_verbose=test_verbose, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                     dgl_shift=test_shift, dgl_tol=tol, dgl_memory=memory, &
                                     dgl_memory_unit=memory_unit, dgl_info=info)
        else
            call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                                     dgl_verbose=test_verbose, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                     dgl_shift=test_shift, dgl_tol=tol, dgl_memory=memory, &
                                     dgl_memory_unit=memory_unit, dgl_info=info)
        end if
!
        same = .false.
        if (ok .and. info .eq. dgl_success) then
            if (generalized) then
                same = compare_eigs(n, n_targ, eig, evec, "Symmetric Generalized")
            else
                same = compare_eigs(n, n_targ, eig, evec, "Symmetric diagonalization")
            end if
            call dump_eigpairs(lutest, n, n_targ, eig, evec, trim(label))
        end if
        call check_result(trim(label), ok, info, same)
!
    end subroutine test_davidson
!
    subroutine test_lobpcg(generalized)
        implicit none
        logical, intent(in) :: generalized
        integer(ip) :: info
        logical :: same
        character(len=120) :: label
!
        if (generalized) then
            label = test_label("Generalized LOBPCG")
        else
            label = test_label("LOBPCG")
        end if
        call init_guess(n, evec)
!
        write (6, f_string) 'testing '//trim(label)//':'
        if (generalized) then
            call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
                                   dgl_verbose=test_verbose, dgl_max_iter=max_iter, dgl_shift=test_shift, &
                                   dgl_tol=tol, dgl_memory=memory, dgl_memory_unit=memory_unit, dgl_info=info, &
                                   dgl_precnd_shift=lobpcg_precnd_shift)
        else
            call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                                   dgl_verbose=test_verbose, dgl_max_iter=max_iter, dgl_shift=test_shift, &
                                   dgl_tol=tol, dgl_memory=memory, dgl_memory_unit=memory_unit, dgl_info=info, &
                                   dgl_precnd_shift=lobpcg_precnd_shift)
        end if
!
        same = .false.
        if (ok .and. info .eq. dgl_success) then
            if (generalized) then
                same = compare_eigs(n, n_targ, eig, evec, "Symmetric Generalized")
            else
                same = compare_eigs(n, n_targ, eig, evec, "Symmetric diagonalization")
            end if
            call dump_eigpairs(lutest, n, n_targ, eig, evec, trim(label))
        end if
        call check_result(trim(label), ok, info, same)
!
    end subroutine test_lobpcg
!
    subroutine test_nonsym_davidson(side, matrix)
!
! test non-symmetric davidson for the given side ("R ", "L " or "LR") and test matrix
!
        implicit none
        character(len=2), intent(in) :: side
        integer(ip), intent(in) :: matrix
        integer(ip) :: info
        logical :: same
        character(len=120) :: label, reference
        procedure(dgl_matvec), pointer :: matvec_r, matvec_l
        procedure(dgl_precnd), pointer :: precnd
!
        select case (matrix)
        case (close_eigs)
            label = "non-symmetric Davidson ("//trim(side)//", nearly degenerate)"
            reference = "Nearly degenerate nonsymmetric diagonalization"
            matvec_r => arx_close
            matvec_l => alx_close
            precnd => dx_close
            if (weak_precnd) precnd => dx_close_weak
        case (cplx_eigs)
            label = "non-symmetric Davidson ("//trim(side)//", complex pair)"
            reference = "Complex pair nonsymmetric diagonalization"
            matvec_r => arx_cplx
            matvec_l => alx_cplx
            precnd => dx
        case default
            label = "non-symmetric Davidson ("//trim(side)//")"
            reference = "Non Symmetric diagonalization"
            matvec_r => arx
            matvec_l => alx
            precnd => dx
        end select
        label = test_label(trim(label))
!
        call init_guess(n, evec)
        call init_guess(n, evec_2)
!
        write (6, f_string) 'testing '//trim(label)//':'
        if (side .eq. "LR") then
            call dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, eig, evec, ok, &
                                           evec_2=evec_2, dgl_verbose=test_verbose, dgl_max_iter=max_iter, &
                                           dgl_dav_iter=dav_iter, dgl_shift=test_shift, dgl_tol=tol, &
                                           dgl_memory=memory, dgl_memory_unit=memory_unit, dgl_info=info)
        else
            call dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, eig, evec, ok, &
                                           dgl_verbose=test_verbose, dgl_max_iter=max_iter, &
                                           dgl_dav_iter=dav_iter, dgl_shift=test_shift, dgl_tol=tol, &
                                           dgl_memory=memory, dgl_memory_unit=memory_unit, dgl_info=info)
        end if
!
! in evec: right eigenvectors for "R " and "LR", left ones for "L ". in evec_2: left ones for "LR"
!
        same = .false.
        if (ok .and. info .eq. dgl_success) then
            select case (side)
            case ("R ")
                same = compare_eigs(n, n_targ, eig, evec, trim(reference)//", Right")
            case ("L ")
                same = compare_eigs(n, n_targ, eig, evec, trim(reference)//", Left")
            case ("LR")
                same = compare_eigs(n, n_targ, eig, evec, trim(reference)//", Right")
                same = compare_eigs(n, n_targ, eig, evec_2, trim(reference)//", Left") .and. same
                same = check_biortho(n, n_targ, evec_2, evec) .and. same
            end select
            call dump_eigpairs(lutest, n, n_targ, eig, evec, trim(label))
        end if
        call check_result(trim(label), ok, info, same)
!
    end subroutine test_nonsym_davidson
!
    subroutine test_smogd()
        implicit none
        integer(ip) :: info
        logical :: same
        real(dp), allocatable :: eig_lr(:), evec_lr(:, :)
        character(len=120) :: label
!
        allocate (evec_lr(2*n, n_max), eig_lr(n_max))
        label = test_label("SMOGD")
        call init_guess(2*n, evec_lr)
!
        write (6, f_string) 'testing '//trim(label)//':'
        call dgl_smogd_driver(2*n, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, eig_lr, evec_lr, ok, &
                              dgl_verbose=test_verbose, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                              dgl_tol=tol, dgl_memory=memory, dgl_memory_unit=memory_unit, dgl_info=info)
        same = .false.
        if (ok .and. info .eq. dgl_success) then
            same = compare_eigs(2*n, n_targ, eig_lr, evec_lr, "Linear response")
            call dump_eigpairs(lutest, n*2, n_targ, eig_lr, evec_lr, trim(label))
        end if
        call check_result(trim(label), ok, info, same)
!
        deallocate (evec_lr, eig_lr)
!
    end subroutine test_smogd
!
    subroutine test_dav_iter()
!
! a value of dav_iter smaller than the minimum is replaced by the minimum; a value such
! that the expansion space would be larger than the problem is reduced.
!
        implicit none
        integer(ip) :: info, i_dav
        integer(ip), parameter :: dav_iters(2) = [1_ip, 200_ip]
        logical :: same
        real(dp), allocatable :: eig_lr(:), evec_lr(:, :)
        character(len=16) :: dav_label
!
        guess_kind = simple_guess
        allocate (evec_lr(2*n, n_max), eig_lr(n_max))
        do i_dav = 1, size(dav_iters)
            write (dav_label, "('dav_iter = ', i0)") dav_iters(i_dav)
!
            call init_guess(n, evec)
            call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, dgl_dav_iter=dav_iters(i_dav), &
                                     dgl_tol=tol, dgl_max_iter=max_iter, dgl_verbose=.true., dgl_info=info)
            same = .false.
            if (ok .and. info .eq. dgl_success) same = compare_eigs(n, n_targ, eig, evec, "Symmetric diagonalization")
            call check_result("Davidson, "//trim(dav_label), ok, info, same)
!
            call init_guess(n, evec)
            call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "R ", eig, evec, ok, &
                                           dgl_dav_iter=dav_iters(i_dav), dgl_tol=tol, dgl_max_iter=max_iter, &
                                           dgl_verbose=.true., dgl_info=info)
            same = .false.
            if (ok .and. info .eq. dgl_success) &
                same = compare_eigs(n, n_targ, eig, evec, "Non Symmetric diagonalization, Right")
            call check_result("non-symmetric Davidson, "//trim(dav_label), ok, info, same)
!
            call init_guess(2*n, evec_lr)
            call dgl_smogd_driver(2*n, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, eig_lr, evec_lr, ok, &
                                  dgl_dav_iter=dav_iters(i_dav), dgl_tol=tol, dgl_max_iter=max_iter, &
                                  dgl_verbose=.true., dgl_info=info)
            same = .false.
            if (ok .and. info .eq. dgl_success) same = compare_eigs(2*n, n_targ, eig_lr, evec_lr, "Linear response")
            call check_result("SMOGD, "//trim(dav_label), ok, info, same)
        end do
        deallocate (evec_lr, eig_lr)
!
    end subroutine test_dav_iter
!
    subroutine test_input_errors()
!
! check that invalid input is reported through dgl_info, without stopping the program
!
        implicit none
        integer(ip) :: info
        procedure(), pointer :: null_p => null()
!
        write (6, f_string) 'testing the handling of invalid input:'
!
        call init_eigenpairs(n, n_max, eig, evec, ok, .false.)
        call dgl_davidson_driver(n, 0_ip, n_max, ax, dx, eig, evec, ok, dgl_info=info)
        call check_error("n_targ = 0", info, dgl_err_input)
        call dgl_davidson_driver(n, n_max + 1, n_max, ax, dx, eig, evec, ok, dgl_info=info)
        call check_error("n_targ > n_max", info, dgl_err_input)
        call dgl_davidson_driver(2*n_max - 1, n_targ, n_max, ax, dx, eig, evec, ok, dgl_info=info)
        call check_error("too many eigenvalues requested", info, dgl_err_input)
        call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, dgl_tol=0.0_dp, dgl_info=info)
        call check_error("tol = 0", info, dgl_err_input)
        call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=null_p, dgl_info=info)
        call check_error("non associated metvec (Davidson)", info, dgl_err_input)
        call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, dgl_max_iter=0_ip, dgl_info=info)
        call check_error("max_iter = 0", info, dgl_err_input)
        call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=null_p, dgl_info=info)
        call check_error("non associated metvec (LOBPCG)", info, dgl_err_input)
        call dgl_lobpcg_driver(3*n_max - 1, n_targ, n_max, ax, dx, eig, evec, ok, dgl_info=info)
        call check_error("too many eigenvalues requested (LOBPCG)", info, dgl_err_input)
        call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "XX", eig, evec, ok, dgl_info=info)
        call check_error("invalid side", info, dgl_err_input)
        call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "LR", eig, evec, ok, dgl_info=info)
        call check_error("side = LR without evec_2", info, dgl_err_input)
        call init_eigenpairs(n, n_max, eig, evec, ok, .false.)
        call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "R ", eig, evec, ok, evec_2=evec_2, &
                                       dgl_tol=tol, dgl_max_iter=max_iter, dgl_info=info)
        call check_error("side = R with an unneeded evec_2 (warning only)", info, dgl_success)
        call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "R ", eig, evec, ok, &
                                       dgl_memory=0_ip, dgl_info=info)
        call check_error("memory = 0", info, dgl_err_input)
        call dgl_smogd_driver(n - 1, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, eig, evec, ok, dgl_info=info)
        call check_error("odd size for SMOGD", info, dgl_err_input)
        call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                                 dgl_memory=1_ip, dgl_memory_unit="KB", dgl_info=info)
        call check_error("not enough memory", info, dgl_err_memory)
        call init_eigenpairs(n, n_max, eig, evec, ok, .false.)
        call dgl_davidson_driver(n, n_targ, n_max, ax_nan, dx, eig, evec, ok, dgl_info=info)
        call check_error("NaN from matvec (Davidson)", info, dgl_err_input)
        call init_eigenpairs(n, n_max, eig, evec, ok, .false.)
        call dgl_lobpcg_driver(n, n_targ, n_max, ax_nan, dx, eig, evec, ok, dgl_info=info)
        call check_error("NaN from matvec (LOBPCG)", info, dgl_err_input)
        call init_eigenpairs(n, n_max, eig, evec, ok, .false.)
        call dgl_davidson_nosym_driver(n, n_targ, n_max, ax_nan, ax_nan, dx, "R ", eig, evec, ok, dgl_info=info)
        call check_error("NaN from matvec (non-symmetric Davidson)", info, dgl_err_input)
        call dgl_smogd_driver(n, n_targ, n_max, ax_nan, ax_nan, spdx, smdx, lrprc, eig, evec, ok, dgl_info=info)
        call check_error("NaN from matvec (SMOGD)", info, dgl_err_input)
!
    end subroutine test_input_errors
!
    subroutine test_dependent_residuals()
!
! the new vectors cannot all be orthonormalized: the dependent ones have to be replaced
!
        implicit none
        integer(ip), parameter :: n_t = 200, n_targ_t = 2, n_max_t = 4
        real(dp), parameter :: ref_std(n_targ_t) = [9.9004942533754781e-01_dp, 1.9999506574411643e+00_dp]
        real(dp), parameter :: ref_gen(n_targ_t) = [9.9749058444760297e-01_dp, 1.9999999999999993e+00_dp]
        real(dp) :: eig_t(n_max_t), evec_t(n_t, n_max_t), evec_t2(n_t, n_max_t)
        procedure(dgl_matvec), pointer :: mx_tri_p
        integer(ip) :: info, i_test, k
        logical :: generalized
        character(len=80) :: label
!
        mx_tri_p => mx_tri
        do i_test = 1, 6
            evec_t = zero
            do k = 1, n_max_t
                evec_t(k, k) = one
            end do
            evec_t2 = evec_t
            generalized = i_test .eq. 2 .or. i_test .eq. 4
            select case (i_test)
            case (1)
                label = "Davidson"
                call dgl_davidson_driver(n_t, n_targ_t, n_max_t, ax_tri, dx_tri, eig_t, evec_t, ok, &
                                         dgl_tol=1.0e-9_dp, dgl_info=info)
            case (2)
                label = "Generalized Davidson"
                call dgl_davidson_driver(n_t, n_targ_t, n_max_t, ax_tri, dx_tri, eig_t, evec_t, ok, &
                                         metvec=mx_tri_p, dgl_tol=1.0e-9_dp, dgl_info=info)
            case (3)
                label = "LOBPCG"
                call dgl_lobpcg_driver(n_t, n_targ_t, n_max_t, ax_tri, dx_tri, eig_t, evec_t, ok, &
                                       dgl_tol=1.0e-9_dp, dgl_info=info)
            case (4)
                label = "Generalized LOBPCG"
                call dgl_lobpcg_driver(n_t, n_targ_t, n_max_t, ax_tri, dx_tri, eig_t, evec_t, ok, &
                                       metvec=mx_tri_p, dgl_tol=1.0e-9_dp, dgl_info=info)
            case (5)
                label = "non-symmetric Davidson (R)"
                call dgl_davidson_nosym_driver(n_t, n_targ_t, n_max_t, ax_tri, ax_tri, dx_tri, "R ", eig_t, &
                                               evec_t, ok, dgl_tol=1.0e-9_dp, dgl_info=info)
            case (6)
                label = "non-symmetric Davidson (LR)"
                call dgl_davidson_nosym_driver(n_t, n_targ_t, n_max_t, ax_tri, ax_tri, dx_tri, "LR", eig_t, &
                                               evec_t, ok, evec_2=evec_t2, dgl_tol=1.0e-9_dp, dgl_info=info)
            end select
            label = trim(label)//", linearly dependent residuals"
!
! with a metric, the expansion space can degenerate to the point that the new vectors cannot be
! made orthogonal to it at all (which happens, for this extreme problem, with some BLAS
! libraries): reporting dgl_err_ortho, instead of returning wrong eigenpairs, is then the
! expected behaviour
!
            if (generalized .and. info .eq. dgl_err_ortho) then
                write (6, f_string) trim(label)//": orthogonalization reported as impossible (accepted)"
                n_tests = n_tests + 1
            else if (generalized) then
                call check_result(trim(label), ok, info, maxval(abs(eig_t(:n_targ_t) - ref_gen)) .lt. 1.0e-9_dp)
            else
                call check_result(trim(label), ok, info, maxval(abs(eig_t(:n_targ_t) - ref_std)) .lt. 1.0e-9_dp)
            end if
        end do
    end subroutine test_dependent_residuals
!
    subroutine test_rotated()
!
! the eigenvalues of ax_rot are d_i = i + 0.5 sin(i)^2, which are increasing
!
        implicit none
        real(dp) :: ref(n_targ), res(n)
        integer(ip) :: info, k, i_test
        real(dp) :: err, max_res
        character(len=80) :: label
!
        do k = 1, n_targ
            ref(k) = rot_eigenvalue(k)
        end do
        do i_test = 1, 3
            evec = zero
            select case (i_test)
            case (1)
                label = "Davidson, rotated matrix, shifted preconditioner"
                call dgl_davidson_driver(n, n_targ, n_max, ax_rot, dx_rot, eig, evec, ok, dgl_tol=tol, &
                                         dgl_max_iter=max_iter, dgl_info=info)
            case (2)
                label = "Davidson, rotated matrix, preconditioner without shift"
                call dgl_davidson_driver(n, n_targ, n_max, ax_rot, dx_rot, eig, evec, ok, dgl_tol=tol, &
                                         dgl_max_iter=max_iter, dgl_info=info, dgl_precnd_shift=.false.)
            case (3)
                label = "LOBPCG, rotated matrix, preconditioner without shift (default)"
                call dgl_lobpcg_driver(n, n_targ, n_max, ax_rot, dx_rot, eig, evec, ok, dgl_tol=tol, &
                                       dgl_max_iter=3*max_iter, dgl_info=info)
            end select
            err = maxval(abs(eig(:n_targ) - ref))
            max_res = zero
            do k = 1, n_targ
                call ax_rot(n, 1_ip, evec(:, k), res)
                max_res = max(max_res, norm2(res - eig(k)*evec(:, k))/norm2(evec(:, k)))
            end do
            write (6, "(t3,a,es9.2,a,es9.2)") "max eigenvalue error: ", err, "   max residual: ", max_res
            call check_result(trim(label), ok, info, err .lt. 1.0e-9_dp .and. max_res .lt. 1.0e-6_dp)
        end do
    end subroutine test_rotated
!
    subroutine test_nested()
!
! generalized Davidson with a matvec that calls DiagLib itself
!
        implicit none
        integer(ip) :: info
        logical :: same
!
        write (6, f_string) 'testing Generalized Davidson, DiagLib called in the matvec:'
        call init_eigenpairs(n, n_max, eig, evec, ok, .false.)
        call dgl_davidson_driver(n, n_targ, n_max, ax_nested, dx, eig, evec, ok, metvec=mx_p, &
                                 dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, dgl_tol=tol, &
                                 dgl_memory=memory, dgl_memory_unit=memory_unit, dgl_info=info)
        same = .false.
        if (ok .and. info .eq. dgl_success) same = compare_eigs(n, n_targ, eig, evec, "Symmetric Generalized")
        write (6, "(t3,a,i0,a,i0)") "nested calls: ", nested_calls, ", failed: ", nested_failed
        call check_result("Generalized Davidson, DiagLib called in the matvec", ok, info, &
                          same .and. nested_calls .gt. 0 .and. nested_failed .eq. 0)
    end subroutine test_nested
!
    logical function check_biortho(ld, m, v_l, v_r) result(passed)
!
! left and right eigenvectors have to be biorthonormal, with normalized right eigenvectors
!
        implicit none
        integer(ip), intent(in) :: ld, m
        real(dp), intent(in) :: v_l(ld, m), v_r(ld, m)
        real(dp) :: err_bi, err_norm
        integer(ip) :: i, j
!
        err_bi = zero
        err_norm = zero
        do j = 1, m
            do i = 1, m
                err_bi = max(err_bi, abs(dot_product(v_l(:, i), v_r(:, j)) - merge(one, zero, i .eq. j)))
            end do
            err_norm = max(err_norm, abs(norm2(v_r(:, j)) - one))
        end do
        write (6, "(t3,a,es9.2,a,es9.2)") "max biorthonormality error: ", err_bi, &
            "   max norm error of the right vectors: ", err_norm
        passed = err_bi .lt. 1.0e-10_dp .and. err_norm .lt. 1.0e-10_dp
    end function check_biortho
!
    subroutine check_error(label, info, expected)
        implicit none
        character(len=*), intent(in) :: label
        integer(ip), intent(in) :: info, expected
!
        n_tests = n_tests + 1
        if (info .eq. expected) then
            write (6, "('-- ', a, ': PASSED (info = ', i0, ')', /)") label, info
        else
            n_failed = n_failed + 1
            write (6, "('-- ', a, ': FAILED (info = ', i0, ', expected ', i0, ')', /)") label, info, expected
        end if
    end subroutine check_error
!
end program test_fortran
