submodule(dgl_interface) dgl_smogd
    use dgl_global_utils
    use dgl_orthogonalizations, only: b_ortho, b_ortho_vs_x
    use dgl_minor_utils
!
    implicit none
!
contains
!
    module subroutine dgl_smogd_driver(n2, n_targ, n_max, apbmul, ambmul, &
                                       spdmul, smdmul, lrprec, eig, evec, ok, &
                                       dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                       dgl_memory, dgl_memory_unit, dgl_info)
!
! the arguments are documented in the interface, in dgl_interface.f90. They are repeated
! here, instead of using "module procedure", so that all compilers check the calls to
! the user-supplied routines.
!
        implicit none
        integer(dgl_int), intent(in) :: n2
        integer(dgl_int), intent(in) :: n_targ
        integer(dgl_int), intent(in) :: n_max
        real(dgl_real), dimension(n_max), intent(inout) :: eig
        real(dgl_real), dimension(n2, n_max), intent(inout) :: evec
        logical, intent(inout) :: ok
        procedure(dgl_matvec) :: apbmul
        procedure(dgl_matvec) :: ambmul
        procedure(dgl_matvec) :: spdmul
        procedure(dgl_matvec) :: smdmul
        procedure(dgl_smogd_precnd) :: lrprec
        logical, optional, intent(in) :: dgl_verbose
        integer(dgl_int), optional, intent(in) :: dgl_dav_iter
        integer(dgl_int), optional, intent(in) :: dgl_max_iter
        integer(dgl_int), optional, intent(in) :: dgl_memory
        character(len=2), optional, intent(in) :: dgl_memory_unit
        real(dgl_real), optional, intent(in) :: dgl_tol
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
        real(dp) :: tol
        character(len=2) :: memory_unit
!
! dimension of the halved problem, the one we are actually solving
!
        integer(ip) :: n
!
! actual expansion space size and total dimension
!
        integer(ip) :: lda, lda2
!
! number of active vectors at a given iteration, and indices to access them
!
        integer(ip) :: n_act, ind, i_beg
!
! current size and total dimension of the expansion space
!
        integer(ip) :: ld_current
!
! number of large arrays that will be allocated
!
        integer(ip) :: n_arrs
!
! number of frozen (i.e. converged) vectors
!
        integer(ip) :: n_frozen
!
        integer(ip) :: it, i_eig
!
        real(dp) :: sqrtn, tol_rms, tol_max
!
! arrays to control convergence and orthogonalization
!
        logical, allocatable :: done(:)
!
! expansion spaces, residuals and their norms.
!
        real(dp), allocatable :: vp(:, :), vm(:, :), lvp(:, :), lvm(:, :), bvp(:, :), bvm(:, :)
        real(dp), allocatable :: rp(:, :), rm(:, :), r_norm(:, :)
!
! eigenvectors of the reduced problem and components of the ritz vectors:
!
        real(dp), allocatable :: up(:, :), um(:, :), eigp(:, :), eigm(:, :), bp(:, :), bm(:, :)
!
! subspace matrix and eigenvalues.
!
        real(dp), allocatable :: s_copy(:, :), s_red_2(:, :), e_red(:)
        real(dp), allocatable :: s_red(:, :)
!
! Scratch vector to avoid recomputation of reduced matrix
!
        real(dp), allocatable :: scratch(:, :)
!
! ================
! START EXECUTION
! ================
!
        ok = .false.
!
! Stupidity checks
!
        if (mod(n2, 2_ip) .ne. 0) call dgl_error(ctx,  &
            "Size of the total problem is not even, something is really wrong with your input", dgl_err_input)
!
        if (4*n_max .ge. n2) call dgl_error(ctx,  &
            "Requested more than half of the total number of eigenvalues: expansions space would break down!", &
            dgl_err_input)
!
! Parse optional arguments
!
        verbose_in = .false.; if (present(dgl_verbose)) verbose_in = dgl_verbose
        max_iter = 100; if (present(dgl_max_iter)) max_iter = dgl_max_iter
        dav_iter = 25; if (present(dgl_dav_iter)) dav_iter = dgl_dav_iter
        tol = 1.e-7_dp; if (present(dgl_tol)) tol = dgl_tol
        memory = 80; if (present(dgl_memory)) memory = dgl_memory
        memory_unit = "MB"; if (present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
! check the input
!
        call dgl_check_input(ctx, n2, n_targ, n_max, max_iter, tol, memory)
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
!
! compute the actual size of the expansion space, checking that
! the input makes sense.
! no expansion space smaller than dav_iter = 10 is deemed acceptable.
!
        n = n2/2
        lda = dav_iter*n_max
        lda2 = 2*lda
        if (lda .ge. n) then
            if (verbose_in) call dgl_warning("Expansion space is larger than the dimension of the problem. "// &
                                             "Reducing size to avoid Rouché-Capelli failure")
            lda = n - 1
        end if
!
! compute the number of large vectors that will be allocated
! to later exstimate required memory in dgl_init
!
        n_arrs = lda*6 + n_max*6
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
! matrix-multiplied vectors and the residual:
!
        call mallocate(ctx, n, lda, vp)
        call mallocate(ctx, n, lda, vm)
        call mallocate(ctx, n, lda, lvp)
        call mallocate(ctx, n, lda, lvm)
        call mallocate(ctx, n, lda, bvp)
        call mallocate(ctx, n, lda, bvm)
        call mallocate(ctx, n, n_max, rp)
        call mallocate(ctx, n, n_max, rm)
!
! allocate memory for convergence check
!
        call mallocate(ctx, n_max, done)
        call mallocate(ctx, 2_ip, n_max, r_norm)
!
! allocate memory for the reduced matrix and its eigenvalues:
!
        call mallocate(ctx, lda, lda, s_copy)
        call mallocate(ctx, lda, lda, s_red_2)
        call mallocate(ctx, lda, lda, s_red)
        call mallocate(ctx, lda2, e_red)
!
! allocate memory for the plus and minus eigenvector components:
!
        call mallocate(ctx, lda, n_max, up)
        call mallocate(ctx, lda, n_max, um)
        call mallocate(ctx, n, n_max, eigp)
        call mallocate(ctx, n, n_max, eigm)
        call mallocate(ctx, n, n_max, bp)
        call mallocate(ctx, n, n_max, bm)
!
        call mallocate(ctx, n_max, lda, scratch)
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
        t_diag = zero
        t_ortho = zero
        t_mv = zero
        t_tot = zero
        vp = zero
        vm = zero
        bvp = zero
        bvm = zero
        lvp = zero
        lvm = zero
        ok = .false.
        done = .false.
!
        call get_time(t_tot)
!
! move the guess into the expansion space.
! as for the other drivers, if no guess is provided (evec is zero), a random one is used.
! the same is done for the individual plus and minus components that vanish
! (e.g., if Y = Z), as they cannot be orthonormalized.
!
        if (dnrm2(n2*n_max, evec, 1_ip) .lt. num_thresh) call random_number(evec)
        do i_eig = 1, n_max
            vp(:, i_eig) = evec(1:n, i_eig) + evec(n + 1:n2, i_eig)
            vm(:, i_eig) = evec(1:n, i_eig) - evec(n + 1:n2, i_eig)
            if (dnrm2(n, vp(:, i_eig), 1_ip) .lt. num_thresh) call random_number(vp(:, i_eig))
            if (dnrm2(n, vm(:, i_eig), 1_ip) .lt. num_thresh) call random_number(vm(:, i_eig))
        end do
!
! initialize the counters
!
        n_act = n_max
        ind = 1
        i_beg = 1
        ld_current = 0
!
! main loop:
!
1030    format(t5, 'SMO-GD iterations (tol=', d10.2, '):', /, &
               t5, '------------------------------------------------------------------', /, &
               t7, '  iter  root              eigenvalue', '         rms         max ok', /, &
               t5, '------------------------------------------------------------------')
1040    format(t9, i4, 2x, i4, f24.12, 2d12.4, l3)
!
        if (verbose_in) write (6, 1030) tol
!
        do it = 1, max_iter
!
! update the size of the expansion space.
!
            ld_current = ld_current + n_act
!
! perform this iteration's matrix-vector multiplications:
!
            call get_time(t1)
!
            call apbmul(n, n_act, vp(1, i_beg), lvp(1, i_beg))
            call ambmul(n, n_act, vm(1, i_beg), lvm(1, i_beg))
!
            call get_time(t2)
            t_mv = t_mv + t2 - t1
!
            call get_time(t1)
!
            call b_ortho(ctx, n, n_act, vp(1, i_beg), lvp(1, i_beg))
            call b_ortho(ctx, n, n_act, vm(1, i_beg), lvm(1, i_beg))
!
            call get_time(t2)
            t_ortho = t_ortho + t2 - t1
!
            call get_time(t1)
!
            call spdmul(n, n_act, vp(1, i_beg), bvm(1, i_beg))
            call smdmul(n, n_act, vm(1, i_beg), bvp(1, i_beg))
!
            call get_time(t2)
            t_mv = t_mv + t2 - t1
!
! update the reduced matrix
!
            call dgemm('t', 'n', ld_current, n_act, n, one, vm, n, bvm(1, i_beg), n, zero, s_red(1, i_beg), lda)
            if (it .gt. 1) then
                call dgemm('t', 'n', n_act, i_beg - 1, n, one, vm(1, i_beg), n, bvm, n, zero, scratch, n_max)
                s_red(i_beg:ld_current, 1:i_beg - 1) = scratch(:n_act, :i_beg - 1)
            end if
!
! save s, and assemble s^t s:
!
            s_copy = s_red
            call dgemm('t', 'n', ld_current, ld_current, ld_current, one, s_copy, lda, s_copy, lda, zero, s_red_2, lda)
!
! diagonalize s^t s
!
            call get_time(t1)
            call dsyev('v', 'u', ld_current, s_red_2, lda, e_red, work, lwork, info)
            call get_time(t2)
            t_diag = t_diag + t2 - t1
            if (info .ne. 0) call dgl_error(ctx, "diagonalization of the reduced matrix failed", dgl_err_lapack)
            if (dgl_failed(ctx)) go to 900
!
! extract the eigenvalues and compute the ritz approximation to the
! eigenvectors
!
            do i_eig = 1, n_max
                eig(i_eig) = sqrt(e_red(ld_current - i_eig + 1))
                up(1:ld_current, i_eig) = s_red_2(1:ld_current, ld_current - i_eig + 1)
            end do
!
! compute the u_- eigenvectors:
!
            call dgemm('n', 'n', ld_current, n_max, ld_current, one, s_copy, lda, up, lda, zero, um, lda)
            do i_eig = 1, n_max
                um(1:ld_current, i_eig) = um(1:ld_current, i_eig)/eig(i_eig)
            end do
!
! asemble the symmetric and antysimmetric combinations (Y+Z) and (Y-Z)
!
            call dgemm('n', 'n', n, n_max, ld_current, one, vp, n, up, lda, zero, eigp, n)
            call dgemm('n', 'n', n, n_max, ld_current, one, vm, n, um, lda, zero, eigm, n)
!
! compute the residuals, and their rms and sup norms:
!
            call dgemm('n', 'n', n, n_max, ld_current, one, bvp, n, um, lda, zero, rp, n)
            call dgemm('n', 'n', n, n_max, ld_current, one, bvm, n, up, lda, zero, rm, n)
            call dgemm('n', 'n', n, n_max, ld_current, one, lvp, n, up, lda, zero, bp, n)
            call dgemm('n', 'n', n, n_max, ld_current, one, lvm, n, um, lda, zero, bm, n)
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
                call daxpy(n, -eig(i_eig), bp(:, i_eig), 1_ip, rp(:, i_eig), 1_ip)
                call daxpy(n, -eig(i_eig), bm(:, i_eig), 1_ip, rm(:, i_eig), 1_ip)
                if (i_eig .gt. n_targ) cycle
                r_norm(1, i_eig) = (dnrm2(n, rp(:, i_eig), 1_ip) + dnrm2(n, rm(:, i_eig), 1_ip))/(eig(i_eig)*sqrt(two)*sqrtn)
                r_norm(2, i_eig) = (maxval(abs(rp(:, i_eig))) + maxval(abs(rm(:, i_eig))))/(sqrt(two)*eig(i_eig))
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
                    write (6, 1040) it, i_eig, one/eig(i_eig), r_norm(:, i_eig), done(i_eig)
                end do
                write (6, *)
            end if
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
                if (verbose_in) write (6, '(t7,a)') 'Restarting davidson.'
!
! put current eigenvectors into the first position of the
! expansion space
!
                vp(:, :n_max) = eigp
                vm(:, :n_max) = eigm
!
                lvp(:, :n_max) = bp
                lvm(:, :n_max) = bm
                call b_ortho(ctx, n, n_max, vp, lvp)
                call b_ortho(ctx, n, n_max, vm, lvm)
!
                call dgemm('n', 'n', n, n_max, ld_current, one, bvp, n, um, lda, zero, bp, n)
                call dgemm('n', 'n', n, n_max, ld_current, one, bvm, n, up, lda, zero, bm, n)
                bvp(:, :n_max) = bp
                bvm(:, :n_max) = bm
!
                s_red = zero
                do i_eig = 1, n_max
                    s_red(i_eig, i_eig) = eig(i_eig)
                end do
!
! initialize indexes back to their starting values
!
                ld_current = n_max
                i_beg = n_max + 1
!
            end if
!
! compute the preconditioned residuals using davidson's procedure
! note that this is done with a user-supplied subroutine, that can
! be generalized to experiment with fancy preconditioners that may
! be more effective than the diagonal one, as in the original
! algorithm.
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
            ind = n_max - n_act + 1
            call lrprec(n, n_act, eig(ind), rp(1, ind), rm(1, ind), vp(1, i_beg), vm(1, i_beg))
!
! orthogonalize the new vectors to the existing ones and then
! orthonormalize them.
!
            call get_time(t1)
!
            call b_ortho_vs_x(ctx, n, ld_current, n_act, vp, lvp, vp(1, i_beg))
            call b_ortho_vs_x(ctx, n, ld_current, n_act, vm, lvm, vm(1, i_beg))
!
            call get_time(t2)
            t_ortho = t_ortho + t2 - t1
            if (dgl_failed(ctx)) go to 900
!
            if (verbose_in) write (6, 1050) n_targ, n_act, n_frozen
!
        end do
!
! assemble the eigenvalues and the eigenvectors (Y, Z) from the latest ritz
! approximation, whether or not the procedure has converged.
!
        if (max_iter .ge. 1) then
            eig = one/eig
            evec(1:n, :) = (eigp + eigm)/two
            evec(n + 1:n2, :) = (eigp - eigm)/two
        end if
!
900     continue
        call get_time(t1)
        t_tot = t1 - t_tot
!
        if (verbose_in) write (6, 1000) t_mv, t_diag, t_ortho, t_tot
!
        call mfree(ctx, work)
        call mfree(ctx, vp)
        call mfree(ctx, vm)
        call mfree(ctx, lvp)
        call mfree(ctx, lvm)
        call mfree(ctx, bvp)
        call mfree(ctx, bvm)
        call mfree(ctx, rp)
        call mfree(ctx, rm)
        call mfree(ctx, r_norm)
        call mfree(ctx, done)
        call mfree(ctx, s_copy)
        call mfree(ctx, s_red_2)
        call mfree(ctx, s_red)
        call mfree(ctx, e_red)
        call mfree(ctx, up)
        call mfree(ctx, um)
        call mfree(ctx, eigp)
        call mfree(ctx, eigm)
        call mfree(ctx, bp)
        call mfree(ctx, bm)
        call mfree(ctx, scratch)
!
        call dgl_check_memleak(ctx)
!
! report the error status (or stop, if dgl_info is not present)
!
999     continue
        if (dgl_failed(ctx)) ok = .false.
        call dgl_return_info(ctx, dgl_info)
!
1000    format(t3, 'timings for SMO-GD (cpu/wall):   ', /, &
               t3, '  matrix-vector multiplications: ', 2f12.4, /, &
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
    end subroutine dgl_smogd_driver

end submodule dgl_smogd
