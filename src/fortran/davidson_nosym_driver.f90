module mod_davidson_nosym_driver
    use dgl_global_utils
    use dgl_orthogonalizations, only: ortho, ortho_vs_x
    use dgl_minor_utils
    use dgl_external_interfaces, only: matvec_, metvec_, precnd_
!
    implicit none
!
contains
!
    subroutine davidson_nosym_driver(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, &
                                     eig, evec_1, ok, evec_2, &
                                     dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                     dgl_shift, dgl_memory, dgl_memory_unit)
!! # Driver for Davidson-Liu non-symmetric diagonalization
!! Non-symmetric davidson diagonalization is commonly encountered in EOM-CC theory.
!! This driver can eveluate both Left and Right eigenvectors.
!! Only standard eigenvalue problems.
!! @note
!! eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
!! @endnote
!!
!! @todo
!! Currently there is a non tested procedure that checks for the possible exchange of eigenvectors
!! throught the iterations. It has to be tested and probably cleaned a little bit.
!! @endtodo
        implicit none
        integer, intent(in) :: n
!! Size of the matrix to be diagonalized
        integer, intent(in) :: n_targ
!! Number of required eigenpairs.
        integer, intent(in) :: n_max
!! Maximum size of the search space. Should be >= n_targ
        character(len=2), intent(in) :: side
!! String to decide which eigenvectors to compute and whether to compute
!! both. Possible values are "R ", "L " or "LR"
        real(dp), dimension(n_max), intent(inout) :: eig
!! Computed eigenvalues
        real(dp), dimension(n, n_max), intent(inout) :: evec_1
!! First set of computed vectors. In input it should contain a guess
!! If side="LR" contains the Right ones
        real(dp), dimension(n, n_max), optional, intent(inout) :: evec_2
!! Second set of computed vectors. First set of converged vectors is used as guess.
!! If side="LR" contains the Left ones
        logical, intent(inout) :: ok
!! True if davidson converged
        procedure(matvec_) :: matvec_r
!! External subroutine that performs the matrix-vector multiplication for
!! right eigenvectors
        procedure(matvec_) :: matvec_l
!! External subroutine that performs the matrix-vector multiplication for
!! left eigenvectors
        procedure(precnd_) :: precnd
!! External subroutine that applies a preconditioner
        logical, optional, intent(in) :: dgl_verbose
!! Verbose mode. Default = .false.
        integer, optional, intent(in) :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
        integer, optional, intent(in) :: dgl_dav_iter
!! Maximum number of iterations before Davidson restart. Default = \(25\)
        integer, optional, intent(in) :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
        character(len=2), optional, intent(in) :: dgl_memory_unit
!! Unit of memory. Default = MB
        real(dp), optional, intent(in) :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
        real(dp), optional, intent(in) :: dgl_shift
!! Diagonal level shifting parameter. Default = \(0.\)
!
! local variables:
! ================
        logical :: verbose_in
        integer :: max_iter, dav_iter, memory
        real(dp) :: tol, shift
        character(len=2) :: memory_unit
!
! expansion space varibles:
! total dimension, current dimension
!
        integer :: lda, ld_current
!
! number of large arrays that will be allocated
!
        integer :: n_arrs
!
! number of active vectors at a given iteration, and indices to access them
!
        integer :: n_act, ind, i_beg
!
! number of frozen (i.e. converged) vectors
!
        integer :: n_frozen
!
! tolerances on residuals norms, used for convergence
!
        real(dp) :: tol_rms, tol_max
!
! iterators and utilities
!
        integer :: it, i_eig
        real(dp) :: sqrtn, tol_im
        real(dp) :: xx(1), yy
        integer :: j
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
        real(dp), allocatable :: copy_evec(:, :)
!
! variables for left, right, or both eigenvectors
!
        character(len=1) :: current_side
        integer :: davidson_runs, current_run
        real(dp) :: eig_r(n_max)
!
! variables for the sorting eigenpairs, since lapack does not.
!
        integer :: max_idx(1)
        logical :: found_im, found_er, double
        logical, allocatable :: mask_overlap(:)
!
        real(dp), allocatable :: overlap(:, :), perm_mat(:, :), evec_temp(:, :), eig_temp(:), overlap_diff(:)
        real(dp) :: overlap_val(n_max, 2), overlap_self(n_max)
        real(dp), allocatable :: perm_temp(:, :)
        integer :: overlap_idx(n_max, 2), k
!
! Scratch vector to avoid recomputation of reduced matrix
!
        real(dp), allocatable :: scratch(:, :)
!
! ================
! START EXECUTION
! ================
!
! Stupidity checks:
!
! Check dimension VS. number of eigs requested
!
        if (n_targ .gt. n_max) call dgl_error( &
            "Number of eigenvalues requested is larger that size of arrays passed")
!
! Parse optional arguments
!! zio pera
        verbose_in = .false.; if (present(dgl_verbose)) verbose_in = dgl_verbose
        max_iter = 100; if (present(dgl_max_iter)) max_iter = dgl_max_iter
        dav_iter = 25; if (present(dgl_dav_iter)) dav_iter = dgl_dav_iter
        tol = 1.e-7_dp; if (present(dgl_tol)) tol = dgl_tol
        shift = 0.e0_dp; if (present(dgl_shift)) shift = dgl_shift
        memory = 80; if (present(dgl_memory)) memory = dgl_memory
        memory_unit = "MB"; if (present(dgl_memory_unit)) memory_unit = dgl_memory_unit
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
                call dgl_error("Missing required array for storing left eigenvectors")
            end if
        case default
            call dgl_error("Invalid value for side, options are: 'L ', 'R ' or 'LR'")
        end select
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
        call dgl_init(n, n_arrs, memory, memory_unit, verbose_in)
!
! start by allocating memory for the various lapack routines
!
        lwork = get_mem_lapack(n, n_max)
        call mallocate(lwork, work)
        call mallocate(lda, tau)
!
! allocate memory for for expansion space, the corresponding
! matrix-multiplied vectors and the residuals
!
        call mallocate(n, lda, space)
        call mallocate(n, lda, aspace)
        call mallocate(n, n_max, residuals)
!
! allocate memory for convergency check
!
        call mallocate(n_max, done)
        call mallocate(2, n_max, r_norm)
!
! allocate memory for the reduced matrix, its eigenvalues with real &
! imaginary parts, and its left & right eigenvectors
!
        call mallocate(lda, lda, a_red)
        call mallocate(lda, lda, a_copy)
        call mallocate(2*lda, e_red_re)
        call mallocate(2*lda, e_red_im)
        call mallocate(lda, lda, evec_red)
        call mallocate(lda, lda, copy_evec)
!
! allocate space for orthogonalization routines
! and mask array for sorting routine
!
        call mallocate(2*n_max, 2*n_max, overlap)
        call mallocate(n_max, overlap_diff)
        call mallocate(2*n_max, 2*n_max, perm_mat)
        call mallocate(2*n_max, 2*n_max, perm_temp)
        call mallocate(lda, n_max, evec_temp)
        call mallocate(lda, eig_temp)
        call mallocate(2*n_max, mask_overlap)
!
        call mallocate(n_max, lda, scratch)
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
        call check_guess(n, n_max, evec_1)
!
! Move guess into the expansion spaces
!
        call dcopy(n*n_max, evec_1, 1, space, 1)
!
! Start Davidson
!
        do current_run = 1, davidson_runs
!
! clean out various quantities
!
            call get_time(t_tot1)
!
            t_tot = zero
            t_diag = zero
            t_ortho = zero
            t_mv = zero
            done = .false.
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
            if (verbose) write (6, 1030) current_side, tol, current_side
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
                if (abs(shift) .gt. num_thresh) call daxpy(n*n_act, shift, space(1, i_beg), 1, aspace(1, i_beg), 1)

!
! get the reduced matrix
!
                select case (current_side)
                case ("R")
                    call dgemm('t', 'n', ld_current, n_act, n, one, space, n, aspace(:, i_beg), n, zero, a_red(1, i_beg), lda)
                    if (it .gt. 1) then
                        call dgemm('t', 'n', n_act, i_beg - 1, n, one, space(:, i_beg), n, aspace, n, zero, scratch, n_max)
                        a_red(i_beg:ld_current, 1:i_beg - 1) = scratch(:n_act, :i_beg - 1)
                    end if
                case ("L")
                    call dgemm('t', 'n', ld_current, n_act, n, one, aspace, n, space(:, i_beg), n, zero, a_red(1, i_beg), lda)
                    if (it .gt. 1) then
                        call dgemm('t', 'n', n_act, i_beg - 1, n, one, aspace(:, i_beg), n, space, n, zero, scratch, n_max)
                        a_red(i_beg:ld_current, 1:i_beg - 1) = scratch(:n_act, :i_beg - 1)
                    end if
                end select
!
                a_copy = a_red
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
                if (info .ne. 0) then
                    call dgl_error("diagonalization of reduced space failed.")
                end if
!
! sort lowest eigenpairs in increasing order in range 2*n_max to ensure that all n_max
! sought eigenpairs are in the range 2*n_max
!
                if (it .eq. 1) then
                    call sort_eigenpairs(ld_current, e_red_re, e_red_im, evec_red, n_max, lda, .true., tol_im)
                else
                    call sort_eigenpairs(ld_current, e_red_re, e_red_im, evec_red, n_max + n_act, lda, .true., tol_im)
                end if
!
! double check for complex contributions in the n_max sought eigenvalues
!
                found_im = .false.
                do j = 1, n_max
                    if (e_red_im(j) .gt. tol_im) found_im = .true.
                end do
!
                if (found_im .and. verbose) then
                    print *
                    call dgl_warning("==========================================")
                    call dgl_warning("Complex contribution in sought eigenvalues")
                    call dgl_warning("==========================================")
                    print *
                end if
!
! compute overlap of old and new eigenvectors in the dimension of the old eigenvectors
! to ensure correct sorting by checking if largest absolute value of column is on the
! diagonal. if not, use the indices of the largest elements to construct a permutation
! matrix to resort the eigenpairs according to the overlap.
!
                if (it .gt. 1) then
!
! compute overlap for the right eigenvectors and extract the index and value of the largest
! and second largest overlap
!
                    call dgemm('t', 'n', 2*n_max, 2*n_max, ld_current, one, copy_evec, lda, evec_red, lda, zero, overlap, 2*n_max)
!
                    found_er = .false.
                    do j = 1, n_max
                        mask_overlap = .true.
                        max_idx = maxloc(abs(overlap(:, j)))
                        overlap_idx(j, 1) = max_idx(1)
                        overlap_self(j) = overlap(j, j)
                        overlap_val(j, 1) = overlap(max_idx(1), j)
                        mask_overlap(max_idx) = .false.
!
! identify if a swapping is necessary
!
                        if (max_idx(1) .ne. j) found_er = .true.
!
! extract index and value of second larges overlap
!
                        max_idx = maxloc(abs(overlap(:, j)), mask=mask_overlap)
                        overlap_idx(j, 2) = max_idx(1)
                        overlap_val(j, 2) = overlap(max_idx(1), j)
                    end do
!
! check if no indices were assigned twice as maximum overlap
!
                    double = .false.
                    do j = 1, n_max
                        do k = 1, n_max
                            if (k .ne. j .and. abs(overlap_idx(j, 1) - overlap_idx(k, 1)) .lt. num_thresh) then
                                double = .true.
                            end if
                        end do
                    end do
!
! try easy fix, by just taking the permutation indexes of the other eigenvector side
!
                    if (double) then
!
! check which second largest overlap is larger and take this indice as max_overlap.
! try for right side only, if no result, dont swap anything and try to continue
! without swapping any eigenvectors
!
                        do j = 1, n_max
                            do k = 1, n_max
                                if (k .ne. j .and. abs(overlap_idx(j, 1) - overlap_idx(k, 1)) .lt. num_thresh) then
                                    if (overlap_val(j, 2) .gt. overlap_val(k, 2)) then
                                        overlap_idx(j, 1) = overlap_idx(j, 2)
                                    else
                                        overlap_idx(k, 1) = overlap_idx(k, 2)
                                    end if
                                end if
                            end do
                        end do
!
! check again if they are the same indices in the max_overlap for the right side.
! if no, then take these indices for the right and left side. if yes, try without
! swapping
!
                        double = .false.
                        do j = 1, n_max
                            do k = 1, n_max
                                if (k .ne. j .and. abs(overlap_idx(j, 1) - overlap_idx(k, 1)) .lt. num_thresh) then
                                    double = .true.
                                end if
                            end do
                        end do
                        if (double) then
                            do j = 1, n_max
                                overlap_idx(j, 1) = j
                            end do
                        end if
                    end if
!
                    if (found_er) then
!
! now permute eigenvectors according to the maxiumum overlap.
! get permutation matrix first
!
                        perm_mat = zero
                        do j = 1, n_max
                            perm_mat(overlap_idx(j, 1), j) = one
                        end do
                        perm_temp = transpose(perm_mat)
!
! now permute left & right eigenvectors and imaginary & real eigenvalues
! note: the use of 't' instead of computing the transpose explicitly obtained in
! a different result
!
                        call dgemm('n', 'n', ld_current, n_max, 2*n_max, one, evec_red, lda, perm_mat, 2*n_max, &
                                   zero, evec_temp, lda)
                        call dcopy(ld_current*n_max, evec_temp, 1, evec_red, 1)
!
                        call dgemv('n', n_max, 2*n_max, one, perm_temp, 2*n_max, e_red_re, 1, zero, eig_temp, 1)
                        call dcopy(n_max, eig_temp, 1, e_red_re, 1)
!
                        call dgemv('n', n_max, 2*n_max, one, perm_temp, 2*n_max, e_red_im, 1, zero, eig_temp, 1)
                        call dcopy(n_max, eig_temp, 1, e_red_im, 1)
!
                    end if
                end if
!
! copy and save the new eigenvectors for the next iteration
!
                copy_evec = evec_red
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
! compute the residuals, and their rms and sup norms
!
                call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, evec_red, lda, zero, residuals, n)
!
                do i_eig = 1, n_targ
!
! if the eigenvalue is already converged, skip it.
!
                    if (done(i_eig)) cycle
!
                    select case (current_run)
                    case (1)
                        call daxpy(n, -eig(i_eig), evec_1(:, i_eig), 1, residuals(:, i_eig), 1)
                    case (2)
                        call daxpy(n, -eig(i_eig), evec_2(:, i_eig), 1, residuals(:, i_eig), 1)
                    end select
                    r_norm(1, i_eig) = dnrm2(n, residuals(:, i_eig), 1)/sqrtn
                    r_norm(2, i_eig) = maxval(abs(residuals(:, i_eig)))
!
                end do
!
! check convergence. lock the first contiguous converged eigenvalues
! by setting the logical array "done" to true
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
! print some information
!
                if (verbose) then
                    do i_eig = 1, n_targ
                        write (6, 1040) it, i_eig, eig(i_eig) - shift, r_norm(:, i_eig), done(i_eig)
                    end do
                    write (6, *)
                end if
!
                if (all(done(1:n_targ))) then
                    ok = .true.
                    exit
                end if
!
! check weather an update is required.
! if not, perform a davidson restart
!
                if (ld_current + n_act .le. lda) then
!
                    i_beg = i_beg + n_act
!
                else
!
                    if (verbose) write (6, '(t7,a)') 'Restarting Davidson'
                    n_act = n_max
!
! put current eigenvectors into the first position of tne
! expansion space. Same thing with their application
!
                    select case (current_run)
                    case (1)
                        call dcopy(n_max*n, evec_1, 1, space, 1)
                        call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, evec_red, lda, zero, evec_1, n)
                        call dcopy(n_max*n, evec_1, 1, aspace, 1)
                    case (2)
                        call dcopy(n_max*n, evec_2, 1, space, 1)
                        call dgemm('n', 'n', n, n_max, ld_current, one, aspace, n, evec_red, lda, zero, evec_2, n)
                        call dcopy(n_max*n, evec_2, 1, aspace, 1)
                    end select
!
! orthogonalize non orthogonal eigenvectors and propagate
! to their application
!
                    call ortho(n, n_max, space, aspace)
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
                call precnd(n, n_act, shift - eig(ind), residuals(1, ind), space(1, i_beg))
!
! orthogonalize the new vectors to the existing ones of the respective other
! space and orthogonalize set of new vectors among each other
!
! Gram-Schmit orthogonalization of residual to the respective subspace
!
                call get_time(t1)
                call ortho_vs_x(n, ld_current, n_act, space, space(1, i_beg), xx, xx)
                call get_time(t2)
!
                t_ortho = t_ortho + t2 - t1
!
                if (verbose) write (6, 1050) n_targ, n_act, n_frozen
!
            end do
!
! end of davidson, print results
!
            call get_time(t_tot2)
            t_tot = t_tot2 - t_tot1
!
! if required, print timings
!
            if (verbose) then
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
                    call dcopy(n*n_max, evec_1, 1, evec_2, 1)
                    call ortho_cd(n, n_max, evec_2, yy, ok)
                    call dcopy(n*n_max, evec_2, 1, space, 1)
!
                    current_side = "L"
!
                else
!
! check if energies are same
!
                    if (maxval(eig_r(:n_targ) - eig(:n_targ)) .gt. tol) then
                        print "(*(d10.2))", eig_r(:n_targ)
                        print "(*(d10.2))", eig(:n_targ)
                        print "(*(d10.2))", eig_r(:n_targ) - eig(:n_targ)
                        call dgl_error("Eigenvalues in the consecutive computation of "// &
                                       "right and left eigenpairs do not match.")
                    end if
!
                end if

            end if
        end do
!
! deallocate memory
!
        call mfree(work)
        call mfree(tau)
        call mfree(space)
        call mfree(aspace)
        call mfree(residuals)
        call mfree(r_norm)
        call mfree(done)
        call mfree(a_red)
        call mfree(a_copy)
        call mfree(e_red_re)
        call mfree(e_red_im)
        call mfree(evec_red)
        call mfree(copy_evec)
        call mfree(overlap)
        call mfree(overlap_diff)
        call mfree(perm_mat)
        call mfree(perm_temp)
        call mfree(evec_temp)
        call mfree(eig_temp)
        call mfree(mask_overlap)
        call mfree(scratch)
!
        call dgl_check_memleak()
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
    end subroutine davidson_nosym_driver

    subroutine sort_eigenpairs(m, w_re, w_im, v, n_want, ldv, ignore, thresh, mask_in)
!
! sort m real & imaginary eigenvalues and right & left eigenvectors of length n
! in decreasing order according to the real eigenvalues in the rang.e of n_want
!
        implicit none
        integer, intent(in) :: m, ldv, n_want
        real(dp), intent(inout) :: w_re(m), w_im(m), v(ldv, m)
        real(dp), intent(in) :: thresh
        logical, intent(in) :: ignore
        logical, optional :: mask_in(m)
!
! local variables
!
        integer :: i, j, idx, min_idx(1), fin
        logical :: mask(m)
!
! define initial mask
!
        if (present(mask_in)) then
            mask = mask_in
        else
            mask = .true.
        end if
!
        do i = 1, n_want
!
! identify minimal value and mask first position for next iteration
!
            min_idx = minloc(w_re, mask=mask)
            idx = min_idx(1)
!
! check complex contribution, if so, move it to with to the last position
! of the array and mask it. search again for lowest eigenvalue and
! continue with that.
!
            if (ignore .and. abs(w_im(idx)) > thresh) then
                fin = m
!
                do j = 1, m
                    if (.not. mask(fin)) then
                        fin = fin - 1
                    else
                        exit
                    end if
                end do
!
                mask(fin) = .false.
!
! do various swaps for double value on last available position fin
!
                call swap_eigenpairs(fin, idx, m, w_re, w_im, v, ldv)
!
! now search again for lowest and find automatically the corresponding
! pair with imaginary contribution
!
                min_idx = minloc(w_re, mask=mask)
                idx = min_idx(1)
            end if
!
            mask(i) = .false.
!
! do various swaps to move minimum value et alii on position i
!
            call swap_eigenpairs(i, idx, m, w_re, w_im, v, ldv)
!
!
        end do
!
    end subroutine sort_eigenpairs

    subroutine swap_eigenpairs(i, j, m, w_re, w_im, v, ldv)
!
! swaps m real & immaginary eigenvalues and eigenvectors of length l of the
! indices i and j with each other
!
        implicit none
        integer, intent(in) :: m, ldv, i, j
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

end module mod_davidson_nosym_driver
