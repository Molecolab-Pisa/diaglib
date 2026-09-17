submodule(dgl_interface) dgl_davidson_nosym
    use dgl_global_utils
    use dgl_orthogonalizations, only: ortho, ortho_vs_x, ortho_gs, biortho_eigvecs
    use dgl_minor_utils
!
    implicit none
!
contains
!
    module subroutine dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, &
                                                eig, evec_1, ok, evec_2, &
                                                dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                                dgl_shift, dgl_memory, dgl_memory_unit, dgl_info, dgl_precnd_shift)
!
! the arguments are documented in the interface, in dgl_interface.f90. They are repeated
! here, instead of using "module procedure", so that all compilers check the calls to
! the user-supplied routines.
!
        implicit none
        integer(dgl_int), intent(in) :: n
        integer(dgl_int), intent(in) :: n_targ
        integer(dgl_int), intent(in) :: n_max
        character(len=2), intent(in) :: side
        real(dgl_real), dimension(n_max), intent(inout) :: eig
        real(dgl_real), dimension(n, n_max), intent(inout) :: evec_1
        real(dgl_real), dimension(n, n_max), optional, intent(inout) :: evec_2
        logical, intent(inout) :: ok
        procedure(dgl_matvec) :: matvec_r
        procedure(dgl_matvec) :: matvec_l
        procedure(dgl_precnd) :: precnd
        logical, optional, intent(in) :: dgl_verbose
        integer(dgl_int), optional, intent(in) :: dgl_max_iter
        integer(dgl_int), optional, intent(in) :: dgl_dav_iter
        integer(dgl_int), optional, intent(in) :: dgl_memory
        character(len=2), optional, intent(in) :: dgl_memory_unit
        real(dgl_real), optional, intent(in) :: dgl_tol
        real(dgl_real), optional, intent(in) :: dgl_shift
        integer(dgl_int), optional, intent(out) :: dgl_info
        logical, optional, intent(in) :: dgl_precnd_shift
!
! local variables:
! ================
        type(dgl_context) :: ctx
!! State of this call: memory bookkeeping, error status and verbosity
        real(dp), allocatable :: work(:)
        integer(ip) :: lwork, info
        real(dp) :: t1(2), t2(2), t_diag(2), t_ortho(2), t_mv(2), t_tot(2)
        logical :: verbose_in, precnd_shift
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
! iterators and utilities
!
        integer(ip) :: it, i_eig
        real(dp) :: sqrtn, tol_im
        real(dp) :: tol_eig
        integer(ip) :: j
!
! convergence of the current run
!
        logical :: run_ok
!
! arrays to control convergence
!
        logical, allocatable :: done(:)
!
! expansion spaces, residuals and their norms
!
        real(dp), allocatable :: space(:, :), aspace(:, :)
        real(dp), allocatable :: residuals(:, :)
        real(dp), allocatable :: r_norm(:, :)
!
! subspace matrix, eigenvalues and real and imaginary parts of the eigenvalues
!
        real(dp), allocatable :: a_red(:, :), a_copy(:, :), e_red_re(:), e_red_im(:)
        real(dp), allocatable :: evec_red(:, :)
!
! variables for left, right, or both eigenvectors
!
        character(len=1) :: current_side
        integer(ip) :: davidson_runs, current_run
        real(dp) :: eig_r(n_max)
!
! variables to sort the ritz pairs and to follow the roots through the iterations:
! convergence of each slot, slots ordered by increasing eigenvalue, scratch
!
        logical :: found_im
        logical, allocatable :: conv(:)
        integer(ip), allocatable :: slot_cand(:), rank_slot(:)
        real(dp), allocatable :: ovl_w(:, :), ovl(:, :), evec_temp(:, :), eig_temp(:)
        integer(ip) :: k, n_act_new
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
! Stupidity checks:
!
! Check dimension VS. number of eigs requested
!
        if (2*n_max .ge. n) call dgl_error(ctx,  &
            "Requested more than half of the total number of eigenvalues: expansions space would break down!", &
            dgl_err_input)
!
! Parse optional arguments
!
        verbose_in = .false.; if (present(dgl_verbose)) verbose_in = dgl_verbose
        precnd_shift = .true.; if (present(dgl_precnd_shift)) precnd_shift = dgl_precnd_shift
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
!
! Check option for problem to solve
!
        select case (side)
        case ("R ")
            current_side = "R"
            davidson_runs = 1
            if (present(evec_2)) then
                call dgl_warning("Present non required second array of vectors")
            end if
        case ("L ")
            current_side = "L"
            davidson_runs = 1
            if (present(evec_2)) then
                call dgl_warning("Present non required second array of vectors")
            end if
        case ("LR")
            current_side = "R"
            davidson_runs = 2
            if (.not. present(evec_2)) then
                call dgl_error(ctx, "Missing required array for storing left eigenvectors", dgl_err_input)
            end if
        case default
            call dgl_error(ctx, "Invalid value for side, options are: 'L ', 'R ' or 'LR'", dgl_err_input)
        end select
        if (dgl_failed(ctx)) go to 999
!
! computing actual size of the expansion space, checking that
! the input makes sense.
!
        lda = dav_iter*n_max
        if (lda .ge. n) then
            if (verbose_in) call dgl_warning("Expansion space is larger than the dimension of the problem. "// &
                                             "Reducing size to avoid Rouché-Capelli failure")
            lda = n - 1
        end if
!
! compute the number of large vectors that will be allocated
! to later exstimate required memory in dgl_init
!
        n_arrs = lda*2 + n_max
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
! allocate memory for for expansion space, the corresponding
! matrix-multiplied vectors and the residuals
!
        call mallocate(ctx, n, lda, space)
        call mallocate(ctx, n, lda, aspace)
        call mallocate(ctx, n, n_max, residuals)
!
! allocate memory for convergency check
!
        call mallocate(ctx, n_max, done)
        call mallocate(ctx, 2_ip, n_max, r_norm)
!
! allocate memory for the reduced matrix, its eigenvalues with real &
! imaginary parts, and its left & right eigenvectors
!
        call mallocate(ctx, lda, lda, a_red)
        call mallocate(ctx, lda, lda, a_copy)
        call mallocate(ctx, 2*lda, e_red_re)
        call mallocate(ctx, 2*lda, e_red_im)
        call mallocate(ctx, lda, lda, evec_red)
!
! allocate space to follow the roots through the iterations
!
        call mallocate(ctx, lda, n_max, ovl_w)
        call mallocate(ctx, n_max, n_max, ovl)
        call mallocate(ctx, lda, n_max, evec_temp)
        call mallocate(ctx, lda, eig_temp)
        call mallocate(ctx, n_max, slot_cand)
        call mallocate(ctx, n_max, rank_slot)
        call mallocate(ctx, n_max, conv)
!
        call mallocate(ctx, n_max, lda, scratch)
        if (dgl_failed(ctx)) go to 900
!
! set the tolerance and compute a useful constant to compute rms norms:
!
        sqrtn = sqrt(real(n, dp))
        tol_rms = tol
        tol_max = 10.0_dp*tol
        tol_im = 1.d-12
!
! check weather we have a guess for the eigenvectors in evec, and
! weather it is orthonormal.
! if evec is zero, create a random guess
!
        call check_guess(ctx, n, n_max, evec_1)
        if (dgl_failed(ctx)) go to 900
!
! Move guess into the expansion spaces
!
        call dcopy(n*n_max, evec_1, 1_ip, space, 1_ip)
!
! Start Davidson. ok is true only if all the runs converge.
!
        ok = .true.
        do current_run = 1, davidson_runs
!
! clean out various quantities
!
            call get_time(t_tot)
!
            t_diag = zero
            t_ortho = zero
            t_mv = zero
            done = .false.
            run_ok = .false.
!
! Re-initialize indexes
!
            n_act = n_max
            i_beg = 1
            ind = 1
!
! initialize the counter for the expansion of the subspace
!
            ld_current = 0
!
! print header
!
            if (verbose_in) write (6, 1030) current_side, tol, current_side
!
! main loop
!
            do it = 1, max_iter
!
! update the size of the expansion space.
!
                ld_current = ld_current + n_act
!
! perform this iteration's matrix-vector multiplications for both
! right and left expansion spaces
!
                call get_time(t1)
                select case (current_side)
                case ("R")
                    call matvec_r(n, n_act, space(1, i_beg), aspace(1, i_beg))
                case ("L")
                    call matvec_l(n, n_act, space(1, i_beg), aspace(1, i_beg))
                end select
                call get_time(t2)
                t_mv = t_mv + t2 - t1

!
! get the reduced matrix
!
                select case (current_side)
                case ("R")
                    call dgemm('t', 'n', ld_current, n_act, n, one, space, n, aspace(1, i_beg), n, zero, a_red(1, i_beg), lda)
                    if (it .gt. 1) then
                        call dgemm('t', 'n', n_act, i_beg - 1, n, one, space(1, i_beg), n, aspace, n, zero, scratch, n_max)
                        a_red(i_beg:ld_current, 1:i_beg - 1) = scratch(:n_act, :i_beg - 1)
                    end if
                case ("L")
                    call dgemm('t', 'n', ld_current, n_act, n, one, aspace, n, space(1, i_beg), n, zero, a_red(1, i_beg), lda)
                    if (it .gt. 1) then
                        call dgemm('t', 'n', n_act, i_beg - 1, n, one, aspace(1, i_beg), n, space, n, zero, scratch, n_max)
                        a_red(i_beg:ld_current, 1:i_beg - 1) = scratch(:n_act, :i_beg - 1)
                    end if
                end select
!
                a_copy = a_red
                call dgl_check_finite(ctx, ld_current, a_copy, lda)
                if (dgl_failed(ctx)) go to 900
!
! diagonalize the reduced matrix
!
                call get_time(t1)
                select case (current_side)
                case ("R")
                    call dgeev('n', 'v', ld_current, a_copy, lda, e_red_re, e_red_im, evec_red, lda, evec_red, lda, &
                               work, lwork, info)
                case ("L")
                    call dgeev('v', 'n', ld_current, a_copy, lda, e_red_re, e_red_im, evec_red, lda, evec_red, lda, &
                               work, lwork, info)
                end select
                call get_time(t2)
                t_diag = t_diag + t2 - t1
                if (info .ne. 0) call dgl_error(ctx, "diagonalization of reduced space failed.", dgl_err_lapack)
                if (dgl_failed(ctx)) go to 900
!
! lapack does not sort the eigenvalues: sort all the ritz pairs by increasing real part,
! moving the ones with a non-negligible imaginary part after the real ones.
!
                call sort_eigenpairs(ld_current, e_red_re, e_red_im, evec_red, ld_current, lda, .true., tol_im)
!
! double check for complex contributions in the n_max sought eigenvalues
!
                found_im = .false.
                do j = 1, n_max
                    if (abs(e_red_im(j)) .gt. tol_im) found_im = .true.
                end do
!
                if (found_im .and. verbose_in) then
                    print *
                    call dgl_warning("==========================================")
                    call dgl_warning("Complex contribution in sought eigenvalues")
                    call dgl_warning("==========================================")
                    print *
                end if
!
! the first n_max ritz pairs are the ones we are looking for. from the second iteration on,
! assign each of them to the slot of the previous ritz vector it overlaps most with, so that
! every slot keeps following the same root even if the order of the eigenvalues changes.
!
                if (it .gt. 1) then
                    select case (current_run)
                    case (1)
                        call follow_roots(n, ld_current, n_max, lda, space, evec_1, evec_red, e_red_re, e_red_im, &
                                          ovl_w, ovl, slot_cand, evec_temp, eig_temp)
                    case (2)
                        call follow_roots(n, ld_current, n_max, lda, space, evec_2, evec_red, e_red_re, e_red_im, &
                                          ovl_w, ovl, slot_cand, evec_temp, eig_temp)
                    end select
                end if
!
! extract the eigenvalues and compute the ritz approximation to the
! eigenvectors
!
                eig = e_red_re(1:n_max)
!
                select case (current_run)
                case (1)
                    call dgemm('n', 'n', n, n_max, ld_current, one, space, n, evec_red, lda, zero, evec_1, n)
                case (2)
                    call dgemm('n', 'n', n, n_max, ld_current, one, space, n, evec_red, lda, zero, evec_2, n)
                end select
!
! compute the residuals, and their rms and sup norms, for all the slots.
! as a root may have moved to a different slot, all the slots are checked at every iteration.
!
                call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, evec_red, lda, zero, residuals, n)
!
                do i_eig = 1, n_max
!
                    select case (current_run)
                    case (1)
                        call daxpy(n, -eig(i_eig), evec_1(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
                    case (2)
                        call daxpy(n, -eig(i_eig), evec_2(:, i_eig), 1_ip, residuals(:, i_eig), 1_ip)
                    end select
                    r_norm(1, i_eig) = dnrm2(n, residuals(:, i_eig), 1_ip)/sqrtn
                    r_norm(2, i_eig) = maxval(abs(residuals(:, i_eig)))
                    conv(i_eig) = r_norm(1, i_eig) .lt. tol_rms .and. &
                                  r_norm(2, i_eig) .lt. tol_max .and. &
                                  it .gt. 1
!
                end do
!
! the n_targ roots we are looking for are the ones with the lowest eigenvalues,
! whatever slot they are in.
!
                call rank_slots(n_max, eig, rank_slot)
!
! lock the first contiguous converged slots by setting the logical array "done" to true
!
                done = .false.
                do i_eig = 1, n_targ
                    if (.not. conv(i_eig)) exit
                    done(i_eig) = .true.
                end do
!
! print some information, with the roots in increasing order
!
                if (verbose_in) then
                    do k = 1, n_targ
                        i_eig = rank_slot(k)
                        write (6, 1040) it, k, eig(i_eig) + shift, r_norm(:, i_eig), conv(i_eig)
                    end do
                    write (6, *)
                end if
!
                if (all(conv(rank_slot(1:n_targ)))) then
                    run_ok = .true.
                    exit
                end if
!
! number of slots that will be expanded in the next iteration
!
                n_act_new = n_max - count(done, kind=ip)
!
! check weather an update is required.
! if not, perform a davidson restart
!
                if (ld_current + n_act_new .le. lda) then
!
                    i_beg = i_beg + n_act
!
                else
!
                    if (verbose_in) write (6, '(t7,a)') 'Restarting Davidson'
                    n_act = n_max
!
! put current eigenvectors into the first position of tne
! expansion space. Same thing with their application.
! evec_1/evec_2 are used as scratch, and then restored.
!
                    select case (current_run)
                    case (1)
                        call dcopy(n_max*n, evec_1, 1_ip, space, 1_ip)
                        call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, evec_red, lda, zero, evec_1, n)
                        call dcopy(n_max*n, evec_1, 1_ip, aspace, 1_ip)
                        call dcopy(n_max*n, space, 1_ip, evec_1, 1_ip)
                    case (2)
                        call dcopy(n_max*n, evec_2, 1_ip, space, 1_ip)
                        call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, evec_red, lda, zero, evec_2, n)
                        call dcopy(n_max*n, evec_2, 1_ip, aspace, 1_ip)
                        call dcopy(n_max*n, space, 1_ip, evec_2, 1_ip)
                    end select
!
! orthogonalize non orthogonal eigenvectors and propagate
! to their application
!
                    call ortho(ctx, n, n_max, space, aspace)
                    if (dgl_failed(ctx)) go to 900
!
! reconstruct first block of the reduced matrix
!
                    a_red = zero
                    select case (current_side)
                    case ("R")
                        call dgemm('t', 'n', n_max, n_max, n, one, space, n, aspace, n, zero, a_red, lda)
                    case ("L")
                        call dgemm('t', 'n', n_max, n_max, n, one, aspace, n, space, n, zero, a_red, lda)
                    end select
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
! the level shift passed to the preconditioner is the lowest eigenvalue among the active slots,
! which, as the roots are followed, is not necessarily the one in the first active slot.
!
                n_act = n_act_new
                n_frozen = n_max - n_act
                ind = n_max - n_act + 1
                call precnd(n, n_act, merge(-minval(eig(ind:n_max)), zero, precnd_shift), residuals(1, ind), &
                            space(1, i_beg))
!
! orthogonalize the new vectors to the existing ones of the respective other
! space and orthogonalize set of new vectors among each other
!
! Gram-Schmit orthogonalization of residual to the respective subspace
!
                call get_time(t1)
                call ortho_vs_x(ctx, n, ld_current, n_act, space, space(1, i_beg))
                call get_time(t2)
!
                t_ortho = t_ortho + t2 - t1
                if (dgl_failed(ctx)) go to 900
!
                if (verbose_in) write (6, 1050) n_targ, n_act, n_frozen
!
            end do
            ok = ok .and. run_ok
!
! return the roots of this run sorted by increasing eigenvalue
! (the residuals are not needed anymore, and are used as scratch)
!
            call rank_slots(n_max, eig, rank_slot)
            eig_temp(1:n_max) = eig(rank_slot)
            eig = eig_temp(1:n_max)
            select case (current_run)
            case (1)
                residuals = evec_1(:, rank_slot)
                evec_1 = residuals
            case (2)
                residuals = evec_2(:, rank_slot)
                evec_2 = residuals
            end select
!
! end of davidson, print results
!
            call get_time(t2)
            t_tot = t2 - t_tot
!
! if required, print timings
!
            if (verbose_in) then
                print *
                write (6, 1100) t_mv, t_diag, t_ortho, t_tot
                print *
                print *
            end if
!
! If only one diagonalization is required we are done,
! otherwise let's use our result as a guess for the left side
!
            if (davidson_runs .eq. 2) then
                if (current_run .eq. 1) then
!
                    eig_r = eig
!
! use evec_1 as guess for evec_2
!
                    call dcopy(n*n_max, evec_1, 1_ip, evec_2, 1_ip)
                    call ortho_gs(ctx, n, n_max, evec_2)
                    if (dgl_failed(ctx)) go to 900
                    call dcopy(n*n_max, evec_2, 1_ip, space, 1_ip)
!
                    current_side = "L"
!
                else
!
! if both runs converged, check that the eigenvalues are the same.
! for a non-symmetric matrix, the error on a ritz value is of the order of the norm of the
! residual, which convergence bounds by sqrt(n)*tol; also allow for round-off errors.
!
                    tol_eig = max(sqrtn*tol, 1.0e3_dp*epsilon(one)*max(one, maxval(abs(eig_r(:n_targ)))))
                    if (ok .and. maxval(abs(eig_r(:n_targ) - eig(:n_targ))) .gt. tol_eig) then
                        print "(*(d10.2))", eig_r(:n_targ)
                        print "(*(d10.2))", eig(:n_targ)
                        print "(*(d10.2))", eig_r(:n_targ) - eig(:n_targ)
                        call dgl_error(ctx, "Eigenvalues in the consecutive computation of "// &
                                       "right and left eigenpairs do not match.", dgl_err_mismatch)
                        go to 900
                    end if
!
! biorthonormalize the converged left and right eigenvectors
!
                    if (ok) then
                        call biortho_eigvecs(ctx, n, n_targ, evec_2, evec_1)
                        if (dgl_failed(ctx)) go to 900
                    end if
!
                end if

            end if
        end do
!
! deallocate memory
!
900     continue
        call mfree(ctx, work)
        call mfree(ctx, space)
        call mfree(ctx, aspace)
        call mfree(ctx, residuals)
        call mfree(ctx, r_norm)
        call mfree(ctx, done)
        call mfree(ctx, a_red)
        call mfree(ctx, a_copy)
        call mfree(ctx, e_red_re)
        call mfree(ctx, e_red_im)
        call mfree(ctx, evec_red)
        call mfree(ctx, ovl_w)
        call mfree(ctx, ovl)
        call mfree(ctx, evec_temp)
        call mfree(ctx, eig_temp)
        call mfree(ctx, slot_cand)
        call mfree(ctx, rank_slot)
        call mfree(ctx, conv)
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
1100    format(t3, '  timings for non-symmetric Davidson (cpu/wall) : ', /, &
               t3, '  matrix-vector multiplications   : ', 2f12.4, /, &
               t3, '  diagonalization                 : ', 2f12.4, /, &
               t3, '  orthogonalization               : ', 2f12.4, /, &
               t3, '                                    ', 24('='), /, &
               t3, '  total                           : ', 2f12.4)
!
1030    format(t5, 'Non-symmetric Davidson iterations for ', a1, '-Vectors (tol=', d10.2, '):', /, &
               t5, '-----------------------------------------------------------------------------', /, &
               t7, '  iter  root          eigenvalue        rms(', a1, ')           max   ok', /, &
               t5, '-----------------------------------------------------------------------------')
!
1040    format(t7, 2i6, f20.15, 2d14.4, l5)
!
1050    format(t5, '----------------------------------------', /, &
               t7, '# target vectors:    ', i4, /, &
               t7, '# new vectors added: ', i4, /, &
               t7, '# converged vectors: ', i4, /, &
               t5, '----------------------------------------')
!
!
    end subroutine dgl_davidson_nosym_driver

    subroutine follow_roots(n, ld, n_max, lda, space, x_old, evec_red, w_re, w_im, &
                            ovl_w, ovl, slot_cand, evec_temp, w_temp)
!
! assign the n_max ritz pairs in the first columns of evec_red (eigenvectors of the reduced
! matrix, expressed in the basis space(:, 1:ld)) to the slots of the ritz vectors of the
! previous iteration, x_old, by maximum overlap, and permute them accordingly.
!
! the overlaps are computed between the actual vectors, i.e.,
!
!   <x_old(:, i) | space evec_red(:, p)> = (space^T x_old)(:, i)^T evec_red(:, p),
!
! so that they are meaningful also after a restart of the expansion space.
! the assignment is greedy: the pairs (slot, ritz vector) are matched in order of
! decreasing absolute (normalized) overlap.
!
        implicit none
        integer(ip), intent(in) :: n, ld, n_max, lda
        real(dp), intent(in) :: space(n, ld), x_old(n, n_max)
        real(dp), intent(inout) :: evec_red(lda, n_max), w_re(n_max), w_im(n_max)
        real(dp), intent(inout) :: ovl_w(lda, n_max), ovl(n_max, n_max), evec_temp(lda, n_max), w_temp(n_max)
        integer(ip), intent(inout) :: slot_cand(n_max)
!
        integer(ip) :: i, p, k, best_i, best_p
        real(dp) :: best, fac
        logical :: used(n_max)
!
! overlaps between the old and the new ritz vectors, normalized
!
        call dgemm('t', 'n', ld, n_max, n, one, space, n, x_old, n, zero, ovl_w, lda)
        call dgemm('t', 'n', n_max, n_max, ld, one, evec_red, lda, ovl_w, lda, zero, ovl, n_max)
        do i = 1, n_max
            do p = 1, n_max
                fac = dnrm2(ld, evec_red(:, p), 1_ip)*dnrm2(n, x_old(:, i), 1_ip)
                if (fac .gt. num_thresh) then
                    ovl(p, i) = abs(ovl(p, i))/fac
                else
                    ovl(p, i) = zero
                end if
            end do
        end do
!
! greedy assignment: slot_cand(i) is the ritz vector assigned to slot i.
! in case of ties, the lowest slot and ritz vector are taken, so that the
! order is not changed when the overlaps do not carry any information.
!
        slot_cand = 0
        used = .false.
        do k = 1, n_max
            best = -one
            best_i = 0
            best_p = 0
            do i = 1, n_max
                if (slot_cand(i) .ne. 0) cycle
                do p = 1, n_max
                    if (used(p)) cycle
                    if (ovl(p, i) .gt. best) then
                        best = ovl(p, i)
                        best_i = i
                        best_p = p
                    end if
                end do
            end do
            slot_cand(best_i) = best_p
            used(best_p) = .true.
        end do
!
! permute the ritz pairs
!
        do i = 1, n_max
            evec_temp(1:ld, i) = evec_red(1:ld, slot_cand(i))
        end do
        evec_red(1:ld, 1:n_max) = evec_temp(1:ld, 1:n_max)
!
        w_temp = w_re(slot_cand)
        w_re = w_temp
        w_temp = w_im(slot_cand)
        w_im = w_temp
!
    end subroutine follow_roots

    subroutine rank_slots(n_max, w, rank_slot)
!
! rank_slot(k) is the slot containing the k-th lowest eigenvalue in w
! (stable: slots with the same eigenvalue keep their order)
!
        implicit none
        integer(ip), intent(in) :: n_max
        real(dp), intent(in) :: w(n_max)
        integer(ip), intent(inout) :: rank_slot(n_max)
!
        integer(ip) :: i, j, tmp
!
        do i = 1, n_max
            rank_slot(i) = i
        end do
        do i = 2, n_max
            tmp = rank_slot(i)
            j = i - 1
            do while (j .ge. 1)
                if (w(rank_slot(j)) .le. w(tmp)) exit
                rank_slot(j + 1) = rank_slot(j)
                j = j - 1
            end do
            rank_slot(j + 1) = tmp
        end do
!
    end subroutine rank_slots

    subroutine sort_eigenpairs(m, w_re, w_im, v, n_want, ldv, ignore, thresh)
!
! sort m real & imaginary eigenvalues and the corresponding eigenvectors in
! increasing order of the real part, so that the lowest n_want eigenpairs are
! in the first n_want positions.
! if ignore is true, eigenvalues whose imaginary part is larger than thresh
! are placed after all the real ones.
!
        implicit none
        integer(ip), intent(in) :: m, ldv, n_want
        real(dp), intent(inout) :: w_re(m), w_im(m), v(ldv, m)
        real(dp), intent(in) :: thresh
        logical, intent(in) :: ignore
!
! local variables
!
        integer(ip) :: i, j, idx
!
! selection sort, limited to the first n_want positions
!
        do i = 1, min(n_want, m)
            idx = i
            do j = i + 1, m
                if (comes_before(j, idx)) idx = j
            end do
            if (idx .ne. i) call swap_eigenpairs(i, idx, m, w_re, w_im, v, ldv)
        end do
!
    contains
!
        logical function comes_before(a, b)
!
! true if the eigenpair in position a must be placed before the one in position b
!
            integer(ip), intent(in) :: a, b
            logical :: complex_a, complex_b
!
            complex_a = ignore .and. abs(w_im(a)) .gt. thresh
            complex_b = ignore .and. abs(w_im(b)) .gt. thresh
            if (complex_a .neqv. complex_b) then
                comes_before = complex_b
            else
                comes_before = w_re(a) .lt. w_re(b)
            end if
        end function comes_before
!
    end subroutine sort_eigenpairs

    subroutine swap_eigenpairs(i, j, m, w_re, w_im, v, ldv)
!
! swaps m real & immaginary eigenvalues and eigenvectors of length l of the
! indices i and j with each other
!
        implicit none
        integer(ip), intent(in) :: m, ldv, i, j
        real(dp), intent(inout) :: w_re(m), w_im(m), v(ldv, m)
!
        real(dp) :: w, v_tmp(ldv)
!
        w = w_re(i)
        w_re(i) = w_re(j)
        w_re(j) = w
!
        w = w_im(i)
        w_im(i) = w_im(j)
        w_im(j) = w

        v_tmp = v(:, i)
        v(:, i) = v(:, j)
        v(:, j) = v_tmp
!
    end subroutine swap_eigenpairs

end submodule dgl_davidson_nosym
