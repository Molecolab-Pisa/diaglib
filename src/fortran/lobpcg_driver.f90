submodule(dgl_interface) dgl_lobpcg
    use dgl_global_utils
    use dgl_orthogonalizations, only: ortho_vs_x, b_ortho, b_ortho_vs_x
    use dgl_minor_utils
!
    implicit none
!
contains
!
    module subroutine dgl_lobpcg_driver(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                                        dgl_verbose, dgl_max_iter, dgl_tol, &
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
        integer(ip) :: max_iter, memory
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
! tolerances on residuals norms, used for convergence
!
        real(dp) :: tol_rms, tol_max
!
! number of active vectors at a given iteration, and indices to access them
!
        integer(ip) :: n_act
!
! indexes to access specific parts of the expansion space
!
        integer(ip) :: ind_x, ind_w, ind_p
!
! varible to determine the type of problem
!
        logical :: generalized
!
! iterators and utilities
!
        integer(ip) :: it, i_eig
        real(dp) :: sqrtn
!
! array to control convergence and orthogonalization
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
        real(dp), allocatable :: a_red(:, :), e_red(:)
        real(dp), allocatable :: u_x(:, :), u_p(:, :), x_new(:, :), ax_new(:, :)
        real(dp), allocatable :: bx_new(:, :)
!
! ================
! START EXECUTION
! ================
!
        ok = .false.
!
! Stupidity check
!
        if (3*n_max .ge. n) call dgl_error(ctx,  &
            "Requested more than a third of the total number of eigenvalues: expansions space would break down!", &
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
        tol = 1.e-7_dp; if (present(dgl_tol)) tol = dgl_tol
        shift = 0.e0_dp; if (present(dgl_shift)) shift = dgl_shift
        memory = 80; if (present(dgl_memory)) memory = dgl_memory !80MBs
        memory_unit = "MB"; if (present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
! check the input
!
        call dgl_check_input(ctx, n, n_targ, n_max, max_iter, tol, memory)
        if (dgl_failed(ctx)) go to 999
!
! set size of the expansion space
!
        lda = 3*n_max
!
        if (generalized) then
            n_arrs = lda*3 + n_max*4
        else
            n_arrs = lda*2 + n_max*3
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
        if (generalized) call mallocate(ctx, n, lda, bspace)

!
! allocate memory for the reduced matrix and its eigenvalues:
!
        call mallocate(ctx, lda, lda, a_red)
        call mallocate(ctx, lda, e_red)
!
! allocate memory for temporary copies of x, ax, and bx:
!
        call mallocate(ctx, n, n_max, x_new)
        call mallocate(ctx, n, n_max, ax_new)
        if (generalized) call mallocate(ctx, n, n_max, bx_new)
!
! allocate memory for convergence check
!
        call mallocate(ctx, n_max, done)
        call mallocate(ctx, 2_ip, n_max, r_norm)
        if (dgl_failed(ctx)) go to 900
!
! clean out:
!
        space = zero
        aspace = zero
        a_red = zero
        if (generalized) bspace = zero
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
! if required, compute b*evec and b-orthogonalize the guess
!
        if (generalized) then
            call get_time(t1)
            call metvec(n, n_max, evec, bx_new)
            call get_time(t2)
            t_mv = t_mv + t2 - t1

            call get_time(t1)
            call b_ortho(ctx, n, n_max, evec, bx_new)
            call get_time(t2)
            t_ortho = t_ortho + t2 - t1
        end if
!
! compute the first eigenpairs by diagonalizing the reduced matrix:
!
        call dcopy(n*n_max, evec, 1_ip, space, 1_ip)
        if (generalized) call dcopy(n*n_max, bx_new, 1_ip, bspace, 1_ip)
!
        call get_time(t1)
        call matvec(n, n_max, space, aspace)
        call get_time(t2)
        t_mv = t_mv + t2 - t1
!
        call dgemm('t', 'n', n_max, n_max, n, one, space, n, aspace, n, zero, a_red, lda)
!
        call get_time(t1)
        call dsyev('v', 'l', n_max, a_red, lda, e_red, work, lwork, info)
        call get_time(t2)
        t_diag = t_diag + t2 - t1
        if (info .ne. 0) call dgl_error(ctx, "diagonalization of the reduced matrix failed", dgl_err_lapack)
        if (dgl_failed(ctx)) go to 900
        eig = e_red(1:n_max)
!
! get the ritz vectors:
!
        call dgemm('n', 'n', n, n_max, n_max, one, space, n, a_red, lda, zero, evec, n)
        call dcopy(n*n_max, evec, 1_ip, space, 1_ip)
        call dgemm('n', 'n', n, n_max, n_max, one, aspace, n, a_red, lda, zero, evec, n)
        call dcopy(n*n_max, evec, 1_ip, aspace, 1_ip)
!
! if required, also get b times the ritz vector:
!
        if (generalized) then
            call dgemm('n', 'n', n, n_max, n_max, one, bspace, n, a_red, lda, zero, evec, n)
            call dcopy(n*n_max, evec, 1_ip, bspace, 1_ip)
        end if
!
! do the first iteration explicitly.
! build the residuals:
!
        call dcopy(n*n_max, aspace, 1_ip, residuals, 1_ip)
        if (generalized) then
            do i_eig = 1, n_max
                call daxpy(n, -eig(i_eig), bspace(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
            end do
        else
            do i_eig = 1, n_max
                call daxpy(n, -eig(i_eig), space(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
            end do
        end if
!
! compute the preconditioned residuals:
!
        ind_x = 1
        ind_w = ind_x + n_max
        call precnd(n, n_max, -eig(ind_x), residuals(1, ind_x), space(1, ind_w))
!
! orthogonalize:
!
        call get_time(t1)
        if (generalized) then
            call b_ortho_vs_x(ctx, n, n_max, n_max, space, bspace, space(1, ind_w))
        else
            call ortho_vs_x(ctx, n, n_max, n_max, space, space(1, ind_w))
        end if
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
        if (dgl_failed(ctx)) go to 900
!
! we are now ready to start the main loop.
! initialize a few parameters
!
        tol_rms = tol
        tol_max = ten*tol
        sqrtn = sqrt(real(n, dp))
        ok = .false.
        done = .false.
        n_act = n_max
!
! x_new holds the latest ritz vectors, returned if the procedure does not converge
!
        call dcopy(n*n_max, space, 1_ip, x_new, 1_ip)
!
1010    format(t5, 'LOBPCG iterations (tol=', d10.2, '):')
1020    format(t5, 'Generalized LOBPCG iterations (tol=', d10.2, '):')
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
! perform this iteration's matrix-vector multiplication and b-orthogonalize
! the latest vectors in case of a generalized problem:
!
            if (generalized) then
                call get_time(t1)
                call metvec(n, n_act, space(1, ind_w), bspace(1, ind_w))
                call get_time(t2)
                t_mv = t_mv + t2 - t1

                call get_time(t1)
                call b_ortho(ctx, n, n_act, space(1, ind_w), bspace(1, ind_w))
                call get_time(t2)
                t_ortho = t_ortho + t2 - t1
            end if
!
            call get_time(t1)
            call matvec(n, n_act, space(1, ind_w), aspace(1, ind_w))
            call get_time(t2)
            t_mv = t_mv + t2 - t1
!
! build the reduced matrix and diagonalize it:
!
            ld_current = n_max + 2*n_act
            if (it .eq. 1) ld_current = 2*n_max
            call dgemm('t', 'n', ld_current, ld_current, n, one, space, n, aspace, n, zero, a_red, lda)
!
            call get_time(t1)
            call dsyev('v', 'l', ld_current, a_red, lda, e_red, work, lwork, info)
            call get_time(t2)
            t_diag = t_diag + t2 - t1
!
! if dsyev failed, print an error message and abort (this should not happen)
!
            if (info .ne. 0) call dgl_error(ctx, "diagonalization of the reduced matrix failed", dgl_err_lapack)
            if (dgl_failed(ctx)) go to 900
            eig = e_red(1:n_max)
!
! update x and ax, and, if required, bx:
!
            call dgemm('n', 'n', n, n_max, ld_current, one, space, n, a_red, lda, zero, x_new, n)
            call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, a_red, lda, zero, ax_new, n)
            if (generalized) then
                call dgemm('n', 'n', n, n_max, ld_current, one, bspace, n, a_red, lda, zero, bx_new, n)
            end if
!
! compute the residuals and their rms and sup norms:
!
            call dcopy(n*n_max, ax_new, 1_ip, residuals, 1_ip)
            do i_eig = 1, n_max
!
! if the eigenvalue is already converged, skip it.
!
                if (done(i_eig)) cycle
!
                if (generalized) then
                    call daxpy(n, -eig(i_eig), bx_new(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
                else
                    call daxpy(n, -eig(i_eig), x_new(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
                end if
                r_norm(1, i_eig) = dnrm2(n, residuals(:, i_eig), 1_ip)/sqrtn
                r_norm(2, i_eig) = maxval(abs(residuals(:, i_eig)))
            end do
!
! only lock the first converged eigenvalues/vectors.
!
            do i_eig = 1, n_max
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
! print some information and check for convergence:
!
            if (verbose_in) then
                do i_eig = 1, n_targ
                    write (6, 1040) it, i_eig, eig(i_eig) + shift, r_norm(:, i_eig), done(i_eig)
                end do
                write (6, *)
            end if
            if (all(done(1:n_targ))) then
                call dcopy(n*n_max, x_new, 1_ip, evec, 1_ip)
                ok = .true.
                exit
            end if
!
! compute the number of active eigenvalues.
! converged eigenvalues and eigenvectors will be locked and kept
! for orthogonalization purposes.
!
            n_act = n_max - count(done, kind=ip)
            ind_x = n_max - n_act + 1
            ind_p = ind_x + n_act
            ind_w = ind_p + n_act
!
! compute the new p and ap vectors only for the active eigenvectors.
! this is done by computing the expansion coefficients u_p of x_new
! -x in the basis of (x,p,w), and then by orthogonalizing then to
! the coefficients u_x of x_new.
!
            call mallocate(ctx, ld_current, n_max, u_x)
            call mallocate(ctx, ld_current, n_act, u_p)
            if (dgl_failed(ctx)) go to 900
!
            call get_coeffs(ctx, lda, ld_current, n_max, n_act, a_red, u_x, u_p)
            if (dgl_failed(ctx)) go to 900
!
! p  = space  * u_p
! ap = aspace * u_p
! bp = bspace * u_p
! note that this is numerically safe, as u_p is orthogonal.
!
            call dgemm('n', 'n', n, n_act, ld_current, one, space, n, u_p, ld_current, zero, evec, n)
            call dcopy(n_act*n, evec, 1_ip, space(1, ind_p), 1_ip)
            call dgemm('n', 'n', n, n_act, ld_current, one, aspace, n, u_p, ld_current, zero, evec, n)
            call dcopy(n_act*n, evec, 1_ip, aspace(1, ind_p), 1_ip)
!
            if (generalized) then
                call dgemm('n', 'n', n, n_act, ld_current, one, bspace, n, u_p, ld_current, zero, evec, n)
                call dcopy(n_act*n, evec, 1_ip, bspace(1, ind_p), 1_ip)
            end if
!
            call mfree(ctx, u_x)
            call mfree(ctx, u_p)
!
! now, move x_new and ax_new into space and aspace.
!
            call dcopy(n*n_max, x_new, 1_ip, space, 1_ip)
            call dcopy(n*n_max, ax_new, 1_ip, aspace, 1_ip)
            if (generalized) then
                call dcopy(n*n_max, bx_new, 1_ip, bspace, 1_ip)
            end if
!
! compute the preconditioned residuals w:
!
            call precnd(n, n_act, -eig(1), residuals(1, ind_x), space(1, ind_w))
!
! orthogonalize w against x and p, and then orthonormalize it:
!
            call get_time(t1)
            if (generalized) then
                call b_ortho_vs_x(ctx, n, n_max + n_act, n_act, space, bspace, space(1, ind_w))
            else
                call ortho_vs_x(ctx, n, n_max + n_act, n_act, space, space(1, ind_w))
            end if
            call get_time(t2)
            t_ortho = t_ortho + t2 - t1
            if (dgl_failed(ctx)) go to 900
!
        end do
!
! if not converged, return the latest ritz vectors (evec is used as scratch)
!
        if (.not. ok) call dcopy(n*n_max, x_new, 1_ip, evec, 1_ip)
!
! deallocate memory and return.
!
900     continue
        call mfree(ctx, u_x)
        call mfree(ctx, u_p)
        call mfree(ctx, work)
        call mfree(ctx, space)
        call mfree(ctx, aspace)
        call mfree(ctx, residuals)
        call mfree(ctx, a_red)
        call mfree(ctx, e_red)
        call mfree(ctx, x_new)
        call mfree(ctx, ax_new)
        call mfree(ctx, done)
        call mfree(ctx, r_norm)
        if (generalized) then
            call mfree(ctx, bspace)
            call mfree(ctx, bx_new)
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
1001    format(t3, 'timings for LOBPCG (cpu/wall): ')
1002    format(t3, 'timings for Generalized LOBPCG (cpu/wall): ')
1000    format(t3, '  matrix-vector multiplications: ', 2f12.4, /, &
               t3, '  diagonalization:               ', 2f12.4, /, &
               t3, '  orthogonalization:             ', 2f12.4, /, &
               t3, '                                 ', 24('='), /, &
               t3, '  total:                         ', 2f12.4)
!

    end subroutine dgl_lobpcg_driver
!
    subroutine get_coeffs(ctx, lda, ld_current, n_max, n_act, a_red, u_x, u_p)
        implicit none
        type(dgl_context), intent(inout) :: ctx
!
! given the eigenvetors of the reduced matrix in a_red, extract
! the expansion coefficients for x_new (u_x) and assemble the
! ones for p_new in u_p.
!
! the coefficients u_p are computed as the difference between the
! coefficients for x_new and x_old, and only the columns associated
! with active eigenvectors are considered.
! u_p is then orthogonalized to u_x: this not only guarantees that
! the p_new vectors will be orthogonal to x_new, but also allows one
! to reuse the ax, aw, and ap vectors to compute ap_new, without
! loosing numerical precision.
!
        integer(ip), intent(in) :: lda, ld_current, n_max, n_act
        real(dp), dimension(lda, lda), intent(in) :: a_red
        real(dp), dimension(ld_current, n_max), intent(inout) :: u_x
        real(dp), dimension(ld_current, n_act), intent(inout) :: u_p
!
        integer(ip) :: ind_x, off_x, i_eig
!
        off_x = n_max - n_act
        ind_x = off_x + 1
!
        u_x(1:ld_current, 1:n_max) = a_red(1:ld_current, 1:n_max)
!
! u_p = u_x for the active vectors only
!
        u_p = u_x(:, ind_x:n_max)
!
! remove the coefficients for x from u_p
!
        do i_eig = 1, n_act
            u_p(off_x + i_eig, i_eig) = u_p(off_x + i_eig, i_eig) - one
        end do
!
! orthogonalize:
!
        call ortho_vs_x(ctx, ld_current, n_max, n_act, u_x, u_p)
!
! all done.
!
        return
!
    end subroutine get_coeffs

end submodule dgl_lobpcg

