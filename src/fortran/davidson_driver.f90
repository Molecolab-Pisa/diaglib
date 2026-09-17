submodule(dgl_interface) dgl_davidson
    use dgl_global_utils
    use dgl_orthogonalizations, only: ortho_vs_x, b_ortho, b_ortho_vs_x
    use dgl_minor_utils
!
    implicit none
!
contains
!
    module subroutine dgl_davidson_driver(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                                          dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                          dgl_shift, dgl_memory, dgl_memory_unit, metvec, dgl_info)
!
! the arguments are documented in the interface, in dgl_interface.f90. They are repeated
! here, instead of using "module procedure", so that all compilers check the calls to
! the user-supplied routines.
!
        implicit none
        integer(dgl_int), intent(in) :: n
        integer(dgl_int), intent(in) :: n_targ
        integer(dgl_int), intent(in) :: n_max
        real(dgl_real), dimension(n_max), intent(inout) :: eig
        real(dgl_real), dimension(n, n_max), intent(inout) :: evec
        logical, intent(inout) :: ok
        procedure(dgl_matvec) :: matvec
        procedure(dgl_precnd) :: precnd
        logical, optional, intent(in) :: dgl_verbose
        integer(dgl_int), optional, intent(in) :: dgl_max_iter
        integer(dgl_int), optional, intent(in) :: dgl_dav_iter
        integer(dgl_int), optional, intent(in) :: dgl_memory
        character(len=2), optional, intent(in) :: dgl_memory_unit
        real(dgl_real), optional, intent(in) :: dgl_tol
        real(dgl_real), optional, intent(in) :: dgl_shift
        procedure(dgl_matvec), pointer, optional :: metvec
        integer(dgl_int), optional, intent(out) :: dgl_info
!
! local variables:
! ================
        type(dgl_context) :: ctx
!! State of this call: memory bookkeeping, error status and verbosity
        real(dp), allocatable :: work(:)
        integer(ip) :: lwork, info
        real(dp) :: t1(2), t2(2), t_diag(2), t_ortho(2), t_mv(2), t_tot(2)
        logical :: verbose_in
        integer(ip) :: max_iter, dav_iter, memory
        real(dp) :: tol, shift
        character(len=2) :: memory_unit
!
! expansion space varibles:
! total dimension, current dimension
!
        integer(ip) :: lda, ld_current
!
! number of large arrays that will be allocated
!
        integer(ip) :: n_arrs
!
! number of active vectors at a given iteration, and indices to access them
!
        integer(ip) :: n_act, ind, i_beg
!
! number of frozen (i.e. converged) vectors
!
        integer(ip) :: n_frozen
!
! tolerances on residuals norms, used for convergence
!
        real(dp) :: tol_rms, tol_max
!
! type of problem (standard or generalized)
!
        logical :: generalized
!
! iterators and utilities
!
        integer(ip) :: it, i_eig
        real(dp) :: sqrtn
!
! array to control convergence
!
        logical, allocatable :: done(:)
!
! expansion spaces, residuals and their norms.
!
        real(dp), allocatable :: space(:, :), aspace(:, :), residuals(:, :), r_norm(:, :)
        real(dp), allocatable :: bspace(:, :)
!
! subspace matrix and eigenvalues.
!
        real(dp), allocatable :: a_red(:, :), a_copy(:, :), e_red(:)
        real(dp), allocatable :: b_evec(:, :)
!
! ================
! START EXECUTION
! ================
!
        ok = .false.
!
! Stupidity checks
!
        if (2*n_max .ge. n) call dgl_error(ctx,  &
            "Requested more than half of the total number of eigenvalues: expansions space would break down!", &
            dgl_err_input)
!
! check what problem we are dealing with
!
        generalized = .false.
        if (present(metvec)) then
            if (.not. associated(metvec)) then
                call dgl_error(ctx, "Non associated pointer to metric-vector product routine", dgl_err_input)
            else
                generalized = .true.
            end if
        end if
!
! Parse optional arguments
!
        verbose_in = .false.; if (present(dgl_verbose)) verbose_in = dgl_verbose
        max_iter = 100; if (present(dgl_max_iter)) max_iter = dgl_max_iter
        dav_iter = 25; if (present(dgl_dav_iter)) dav_iter = dgl_dav_iter
        tol = 1.e-7_dp; if (present(dgl_tol)) tol = dgl_tol
        shift = 0.e0_dp; if (present(dgl_shift)) shift = dgl_shift
        memory = 80; if (present(dgl_memory)) memory = dgl_memory
        memory_unit = "MB"; if (present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
! check the input
!
        call dgl_check_input(ctx, n, n_targ, n_max, max_iter, tol, memory)
!
! no expansion space smaller than dgl_min_dav_iter iterations is deemed acceptable
!
        if (dav_iter .lt. dgl_min_dav_iter) then
            if (verbose_in) call dgl_warning("dav_iter is smaller than the minimum allowed value, "// &
                                             "the minimum value is used instead")
            dav_iter = dgl_min_dav_iter
        end if
        if (dgl_failed(ctx)) go to 999
!
! compute the actual size of the expansion space
!
        lda = dav_iter*n_max
        if (lda .ge. n) then
            if (verbose_in) call dgl_warning("Expansion space is larger than the dimension of the problem. "// &
                                             "Reducing size to avoid Rouché-Capelli failure")
            lda = n - 1
        end if
!
! compute the number of large vectors that will be allocated
! to estimate required memory in dgl_init
!
        if (generalized) then
            n_arrs = lda*3 + n_max*2
        else
            n_arrs = lda*2 + n_max
        end if
        call dgl_init(ctx, n, n_arrs, memory, memory_unit, verbose_in)
        t_tot = zero
        t_diag = zero
        t_ortho = zero
        t_mv = zero
!
! start by allocating memory for the various lapack routines
!
        lwork = get_mem_lapack(lda)
        call mallocate(ctx, lwork, work)
!
! allocate memory for the expansion space, the corresponding
! matrix-multiplied vectors and the residuals:
!
        call mallocate(ctx, n, lda, space)
        call mallocate(ctx, n, lda, aspace)
        call mallocate(ctx, n, n_max, residuals)
!
        if (generalized) then
            call mallocate(ctx, n, lda, bspace)
            call mallocate(ctx, n, n_max, b_evec)
        end if
!
! allocate memory for convergence check
!
        call mallocate(ctx, n_max, done)
        call mallocate(ctx, 2_ip, n_max, r_norm)
!
! allocate memory for the reduced matrix and its eigenvalues:
!
        call mallocate(ctx, lda, lda, a_red)
        call mallocate(ctx, lda, lda, a_copy)
        call mallocate(ctx, lda, e_red)
        if (dgl_failed(ctx)) go to 900
!
! set the tolerances and compute a useful constant to compute rms norms:
!
        sqrtn = sqrt(real(n, dp))
        tol_rms = tol
        tol_max = 10.0_dp*tol
!
! clean out various quantities
!
        space = zero
        aspace = zero
        a_red = zero
        if (generalized) bspace = zero
        ok = .false.
        done = .false.
!
        call get_time(t_tot)
!
! check whether we have a guess for the eigenvectors in evec, and
! whether it is orthonormal.
! if evec is zero, create a random guess.
!
        call check_guess(ctx, n, n_max, evec)
        if (dgl_failed(ctx)) go to 900
!
! move the guess into the expansion space.
!
        call dcopy(n*n_max, evec, 1_ip, space, 1_ip)
!
! initialize the number of active vectors and the associated indices.
!
        n_act = n_max
        ind = 1
        i_beg = 1
!
! initialize the counter for the expansion of the subspace
!
        ld_current = 0
!
! main loop:
!
1010    format(t5, 'Davidson-Liu iterations (tol=', d10.2, '):')
1020    format(t5, 'Generalized Davidson-Liu iterations (tol=', d10.2, '):')
1030    format(t5, '------------------------------------------------------------------', /, &
               t7, '  iter  root              eigenvalue', '         rms         max ok', /, &
               t5, '------------------------------------------------------------------')
1040    format(t9, i4, 2x, i4, f24.12, 2d12.4, l3)
!
        if (verbose_in) then
            if (generalized) then
                write (6, 1020) tol
            else
                write (6, 1010) tol
            end if
            write (6, 1030)
        end if
!
        do it = 1, max_iter
!
! update the size of the expansion space.
!
            ld_current = ld_current + n_act
!
! perform this iteration's matrix-vector multiplication and b-orthogonalize
! the latest vectors in case of a generalized problem:
!
            if (generalized) then
                call get_time(t1)
                call metvec(n, n_act, space(1, i_beg), bspace(1, i_beg))
                call get_time(t2)
                t_mv = t_mv + t2 - t1
!
                call get_time(t1)
                call b_ortho(ctx, n, n_act, space(1, i_beg), bspace(1, i_beg))
                call get_time(t2)
                t_ortho = t_ortho + t2 - t1
            end if
!
            call get_time(t1)
            call matvec(n, n_act, space(1, i_beg), aspace(1, i_beg))
            call get_time(t2)
            t_mv = t_mv + t2 - t1
!
! update the reduced matrix
!
            call dgemm('t', 'n', ld_current, n_act, n, one, space, n, aspace(1, i_beg), n, zero, a_red(1, i_beg), lda)
!
! explicitly putting the first block of
! converged eigenvalues in the reduced matrix
!
            a_copy = a_red
!
! diagonalize the reduced matrix
!
            call get_time(t1)
            call dsyev('v', 'u', ld_current, a_copy, lda, e_red, work, lwork, info)
            call get_time(t2)
            t_diag = t_diag + t2 - t1
            if (info .ne. 0) call dgl_error(ctx, "diagonalization of the reduced matrix failed", dgl_err_lapack)
            if (dgl_failed(ctx)) go to 900
!
! extract the eigenvalues and compute the ritz approximation to the
! eigenvectors
!
            eig = e_red(1:n_max)
!
            call dgemm('n', 'n', n, n_max, ld_current, one, space, n, a_copy, lda, zero, evec, n)
!
! compute the residuals, and their rms and sup norms:
!
            call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, a_copy, lda, zero, residuals, n)
            if (generalized) call dgemm('n', 'n', n, n_max, ld_current, one, bspace, n, a_copy, lda, zero, b_evec, n)
!
            do i_eig = 1, n_max
!
! the residuals of all the non-converged ritz pairs are used to expand the space,
! but convergence is only checked for the n_targ lowest ones.
!
                if (i_eig .le. n_targ) then
                    if (done(i_eig)) cycle
                end if
!
                if (generalized) then
                    call daxpy(n, -eig(i_eig), b_evec(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
                else
                    call daxpy(n, -eig(i_eig), evec(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
                end if
                if (i_eig .gt. n_targ) cycle
                r_norm(1, i_eig) = dnrm2(n, residuals(:, i_eig), 1_ip)/sqrtn
                r_norm(2, i_eig) = maxval(abs(residuals(:, i_eig)))
            end do
!
! check convergence. lock the first contiguous converged eigenvalues
! by setting the logical array "done" to true.
!
            do i_eig = 1, n_targ
                if (done(i_eig)) cycle
                done(i_eig) = r_norm(1, i_eig) .lt. tol_rms .and. &
                              r_norm(2, i_eig) .lt. tol_max .and. &
                              it .gt. 1
                if (.not. done(i_eig)) then
                    done(i_eig + 1:n_max) = .false.
                    exit
                end if
            end do
!
! print some information:
!
            if (verbose_in) then
                do i_eig = 1, n_targ
                    write (6, 1040) it, i_eig, eig(i_eig) + shift, r_norm(:, i_eig), done(i_eig)
                end do
                write (6, *)
            end if
!
! exit if everything is converged
!
            if (all(done(1:n_targ))) then
                ok = .true.
                exit
            end if
!
! check whether an update is required.
! if not, perform a davidson restart.
!
            if (ld_current + n_act .le. lda) then
!
                i_beg = i_beg + n_act
!
            else
!
                if (verbose_in) write (6, '(t7,a,/)') 'Restarting davidson'
!
! put current eigenvectors into the first position of the
! expansion space
!
                space(:, n_max + 1:) = zero
                call dcopy(n_max*n, evec, 1_ip, space, 1_ip)
!
                call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, a_copy, lda, zero, evec, n)
                aspace(:, n_max + 1:) = zero
                call dcopy(n_max*n, evec, 1_ip, aspace, 1_ip)
!
! evec has been used as scratch: put the ritz vectors back
!
                call dcopy(n_max*n, space, 1_ip, evec, 1_ip)
!
                a_red = zero
                do i_eig = 1, n_max
                    a_red(i_eig, i_eig) = e_red(i_eig)
                end do
!
                if (generalized) then
                    bspace(:, n_max + 1:) = zero
                    call dcopy(n_max*n, b_evec, 1_ip, bspace, 1_ip)
                    call b_ortho(ctx, n, n_max, space, bspace)
                end if
!
! initialize indexes back to their starting values
!
                ld_current = n_max
                i_beg = n_max + 1
!
            end if
!
! Update number of searched vectors based on converged
! ones: Locking
!
            n_act = n_max
            n_frozen = 0
            do i_eig = 1, n_targ
                if (done(i_eig)) then
                    n_act = n_act - 1
                    n_frozen = n_frozen + 1
                else
                    exit
                end if
            end do
!
! compute the preconditioned residuals using davidson's procedure
! note that this is done with a user-supplied subroutine, that can
! be generalized to experiment with fancy preconditioners that may
! be more effective than the diagonal one, as in the original
! algorithm.
!
            ind = n_max - n_act + 1
            call precnd(n, n_act, -eig(ind), residuals(1, ind), space(1, i_beg))
!
! orthogonalize the new vectors to the existing ones and then
! orthonormalize them.
!
            call get_time(t1)
            if (generalized) then
                call b_ortho_vs_x(ctx, n, ld_current, n_act, space, bspace, space(1, i_beg))
            else
                call ortho_vs_x(ctx, n, ld_current, n_act, space, space(1, i_beg))
            end if
            call get_time(t2)
            t_ortho = t_ortho + t2 - t1
            if (dgl_failed(ctx)) go to 900
!
            if (verbose_in) write (6, 1050) n_targ, n_act, n_frozen
!
        end do
!
! deallocate memory
!
900     continue
        call mfree(ctx, work)
        call mfree(ctx, space)
        call mfree(ctx, aspace)
        call mfree(ctx, bspace)
        call mfree(ctx, residuals)
        call mfree(ctx, done)
        call mfree(ctx, r_norm)
        call mfree(ctx, a_red)
        call mfree(ctx, a_copy)
        call mfree(ctx, e_red)
        if (generalized) then
            call mfree(ctx, bspace)
            call mfree(ctx, b_evec)
        end if
!
        call dgl_check_memleak(ctx)
!
! if required, print timings
!
        call get_time(t2)
        t_tot = t2 - t_tot
        if (verbose_in) then
            if (generalized) then
                write (6, 1002)
            else
                write (6, 1001)
            end if
            write (6, 1000) t_mv, t_diag, t_ortho, t_tot
        end if
!
! report the error status (or stop, if dgl_info is not present)
!
999     continue
        if (dgl_failed(ctx)) ok = .false.
        call dgl_return_info(ctx, dgl_info)
!
1001    format(t3, 'timings for Davidson-Liu (cpu/wall): ')
1002    format(t3, 'timings for Generalized Davidson-Liu (cpu/wall): ')
1000    format(t3, '  matrix-vector multiplications: ', 2f12.4, /, &
               t3, '  diagonalization:               ', 2f12.4, /, &
               t3, '  orthogonalization:             ', 2f12.4, /, &
               t3, '                                 ', 24('='), /, &
               t3, '  total:                         ', 2f12.4)
!
1050    format(t5, '----------------------------------------', /, &
               t7, '# target vectors:    ', i4, /, &
               t7, '# new vectors added: ', i4, /, &
               t7, '# converged vectors: ', i4, /, &
               t5, '----------------------------------------')
!
    end subroutine dgl_davidson_driver

end submodule dgl_davidson
