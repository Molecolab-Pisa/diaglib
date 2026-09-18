module dgl_orthogonalizations
!* Module for all the orthogonalization procedures.
! These are internal routines of the drivers and are not part of the public interface.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use dgl_global_utils
    implicit none
!
    real(dp), parameter, private :: tol_ortho = two*epsilon(one)
!! Convergence thresholds for orthogonalizations
    real(dp), parameter, private :: max_growth_cd = 1.0e10_dp
!! Largest norm of the transformation applied in one iteration of the Cholesky-based
!! orthonormalization: beyond it, the vectors are considered (numerically) linearly dependent,
!! and the Gram-Schmidt fallback is used
    integer(ip), parameter, private :: max_default_it = 3_ip
!! Number of iterations of the default (Cholesky-based) route of orthonormalize_vs_x after
!! which the Gram-Schmidt fallback is used
    real(dp), parameter, private :: keep_ratio = 0.5_dp
!! In gs_fallback, a projection that keeps less than this fraction of the norm of a vector
!! is repeated
!
contains
!
    subroutine ortho(ctx, n, m, u, w)
!! Orthogonalization routine based on QR decomposition.
!! Orthogonalizes \(m\) vectors of length \(n\) contained in \(u\).
!! \[ u^Tu = \textbf{I} \]
!!
!! This is done by computing U = QR and then by solving the upper
!! triangular system U(ortho)R = U.
!! Using this strategy allows to apply the same linear transformation
!! that orthogonalizes U to a second set of vectors that usually contain
!! the product AU, where A is some matrix that has already been applied
!! to U. This is useful when U and AU are built together without explicitly
!! performing the matrix vector multiplication.
        implicit none
        type(dgl_context), intent(inout) :: ctx
!
        integer(ip), intent(in) :: n
!! Length of the input vectors
        integer(ip), intent(in) :: m
!! Number of input vectors
        real(dp), dimension(n, m), intent(inout) :: u
!! Vectors to orthogonalize
        real(dp), dimension(n, m), optional, intent(inout) :: w
!! Second set of vectors to which we may apply the same transformation
!
! local scratch
! =============
!
        real(dp), allocatable :: v(:, :), tau_qr(:), work_qr(:)
        real(dp) :: lwork_query(1)
        integer(ip) :: lwork_qr, info_qr
!
        call mallocate(ctx, n, m, v)
!
! use a local workspace, so that this routine does not depend on
! the global lapack work arrays allocated by the drivers.
!
        call mallocate(ctx, min(n, m), tau_qr)
        if (dgl_failed(ctx)) go to 100
        v = u
        call dgeqrf(n, m, u, n, tau_qr, lwork_query, -1_ip, info_qr)
        lwork_qr = max(1_ip, int(lwork_query(1), ip))
        call mallocate(ctx, lwork_qr, work_qr)
        if (dgl_failed(ctx)) go to 100
        call dgeqrf(n, m, u, n, tau_qr, work_qr, lwork_qr, info_qr)
        if (info_qr .ne. 0) then
            call dgl_error(ctx, "ortho: QR factorization failed", dgl_err_lapack)
            go to 100
        end if
!
        call dtrsm('r', 'u', 'n', 'n', n, m, one, u, n, v, n)
!
        if (present(w)) call dtrsm('r', 'u', 'n', 'n', n, m, one, u, n, w, n)
!
        u = v
!
100     continue
        call mfree(ctx, work_qr)
        call mfree(ctx, tau_qr)
        call mfree(ctx, v)
    end subroutine ortho
!
    subroutine b_ortho(ctx, n, m, u, bu)
!! Subroutine to B-orthogonalize \(m\) vectors of length \(n\)
!! using the Cholesky factorization of their overlap.
!! \[ u^TBu = \textbf{I} \]
!!
!! This is in principle not a good idea, as the \(u^TBu\) matrix can be very
!! ill-conditioned, independently of how bad is \(B\). Also, only works if u is
!! already orthonormal.
!! More details on the orthoganlization scheme are available in the documentation
!! of the [[ortho_cd]] procedure
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip) :: info
!
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of vectors
        real(dp), dimension(n, m), intent(inout) :: u
!! Vectors to B-orthogonalize
        real(dp), dimension(n, m), intent(inout) :: bu
!! Application of an external matrix \(B\) on \(u\)
!
! local variables
! ===============
!
        integer(ip) :: i, j
        real(dp), allocatable :: metric(:, :), sigma(:), u_svd(:, :), vt_svd(:, :), &
                                 temp(:, :), work_svd(:)
        real(dp) :: lwork_query(1)
        integer(ip) :: lwork_svd, info_svd
        real(dp), parameter :: tol_svd = 1.0e-5_dp
        logical, parameter :: use_svd = .false.
!
!
        call mallocate(ctx, m, m, metric)
        if (dgl_failed(ctx)) go to 100
!
        call dgemm('t', 'n', m, m, n, one, u, n, bu, n, zero, metric, m)
        call dgl_check_finite(ctx, m, metric, m)
        if (dgl_failed(ctx)) go to 100
!
        if (use_svd) then
!
! debug option: use svd to b-orthonormalize, by computing
! b**(-1/2)
!
            call mallocate(ctx, m, sigma)
            call mallocate(ctx, m, m, u_svd)
            call mallocate(ctx, m, m, vt_svd)
            call mallocate(ctx, n, m, temp)
            if (dgl_failed(ctx)) go to 100
!
            call dgesvd('a', 'a', m, m, metric, m, sigma, u_svd, m, vt_svd, m, lwork_query, -1_ip, info_svd)
            lwork_svd = max(1_ip, int(lwork_query(1), ip))
            call mallocate(ctx, lwork_svd, work_svd)
            if (dgl_failed(ctx)) go to 100
            call dgesvd('a', 'a', m, m, metric, m, sigma, u_svd, m, vt_svd, m, work_svd, lwork_svd, info_svd)
            call mfree(ctx, work_svd)
!
! compute sigma**(-1/2)
!
            do i = 1, m
                if (sigma(i) .gt. tol_svd) then
                    sigma(i) = 1/sqrt(sigma(i))
                else
                    sigma(i) = zero
                end if
            end do
!
! compute metric ** (-1/2). first, compute sigma ** (-1/2) vt
!
            metric = zero
            do i = 1, m
                do j = 1, m
                    metric(j, i) = metric(j, i) + sigma(j)*vt_svd(j, i)
                end do
            end do
!
! now, multiply for u:
!
            vt_svd = metric
            call dgemm('n', 'n', m, m, m, one, u_svd, m, vt_svd, m, zero, metric, m)
!
! metric contains s ** (-1/2), and projects out directions corresponding
! to pathological singular values.
! orthogonalize u and bu:
!
            call dgemm('n', 'n', n, m, m, one, u, n, metric, m, zero, temp, n)
            u = temp
            call dgemm('n', 'n', n, m, m, one, bu, n, metric, m, zero, temp, n)
            bu = temp
        else
!
! compute the cholesky factorization of the metric.
!
            call dpotrf('l', m, metric, m, info)
            if (info .ne. 0) then
                call dgl_error(ctx, "b_ortho: the metric is not positive definite", dgl_err_ortho)
                go to 100
            end if
!
! get u * l^-T and bu * l^-T
!
            call dtrsm('r', 'l', 't', 'n', n, m, one, metric, m, u, n)
            call dtrsm('r', 'l', 't', 'n', n, m, one, metric, m, bu, n)
        end if
!
100     continue
        call mfree(ctx, sigma)
        call mfree(ctx, u_svd)
        call mfree(ctx, vt_svd)
        call mfree(ctx, temp)
        call mfree(ctx, work_svd)
        call mfree(ctx, metric)
!
    end subroutine b_ortho
!
    subroutine diag_shift(n, shift, a)
!! Add a shift to the diagonal elements of the matrix \(a\)
        implicit none
!
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: shift
        real(dp), dimension(n, n), intent(inout) :: a
!
        integer(ip) :: i
!
        do i = 1, n
            a(i, i) = a(i, i) + shift
        end do
!
        return
    end subroutine diag_shift
!
    subroutine ortho_cd(ctx, n, m, u, growth, ok, max_growth, col_norms)
!! Subroutine to orthogonalize \(m\) vectors of length \(n\)
!! using the Cholesky factorization of their overlap.
!! \[ u^Tu = \textbf{I} \]
!! The metric is computed as \(metric = u^Tu \) and then by computing its cholesky
!! decompositoin \( metric = LL^T \). The orthogonal vectors are obtained then
!! by solving the triangular linear system \( u(ortho)L^T = u \).
!!
!! The metric is scaled by the norms of the vectors, which are the square roots of its diagonal
!! elements, before it is factorized, and the scaling is applied to \(u\) together with the
!! triangular factor: the result does not depend on the norms of the input vectors, and neither
!! does the growth factor, which then only measures their conditioning. The norms are returned
!! in col_norms, if present, so that the caller does not have to compute them.
!!
!! As cholesky decomposition is not the most stable way of orthogonalizing
!! a set of vectors, the orthogonalization is refined iteratively.
!! A conservative estimate of the orthogonalization error is used to
!! assess convergence.
!!
!! This routine returns a growth factor, which can be used in (b_)ortho_vs_x
!! to estimate the orthogonality error introduced by ortho_cd.
!! While it is very unlikely to do so, this routine can fail.
!! The status of this procedure is retured in orther to invoke
!! more robust routines in case of failure.
!! If max_growth is present, the procedure also fails if the norm of the transformation of one
!! iteration would be larger than max_growth, which means that the vectors are (numerically)
!! linearly dependent: that transformation is not applied, so that u still spans the same space.
!! Exactly dependent vectors also make the procedure fail, as the level-shifted iterations, which
!! do not orthonormalize them, are never accepted.
!
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip) :: info
!
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of vectors
        real(dp), dimension(n, m), intent(inout) :: u
!! Vectors to orthogonalize
        real(dp), intent(inout) :: growth
!! Growth factor for the numerical error
        logical, intent(inout) :: ok
!! Status of the procedure in output
        real(dp), optional, intent(in) :: max_growth
!! Maximum growth factor allowed
        real(dp), optional, intent(out) :: col_norms(m)
!! Norms of the vectors in input
!
! local variables
! ===============
!
        integer(ip) :: it, it_micro, i, j
        real(dp) :: error, alpha, shift
        real(dp) :: rcond, l_norm, linv_norm
        logical :: macro_done, micro_done, shifted
        integer(ip), parameter :: maxit = 10
!
! local scratch
! =============
!
        real(dp), allocatable :: metric(:, :), msave(:, :), norms(:)
!
!
! get memory for the metric.
!
        ok = .false.
        call mallocate(ctx, m, m, metric)
        call mallocate(ctx, m, m, msave)
        call mallocate(ctx, m, norms)
        if (dgl_failed(ctx)) go to 100
!
        metric = zero
        macro_done = .false.
!
! assemble the metric
!
        it = 0
        growth = one
        do while (.not. macro_done)
            it = it + 1
!
! ortho_cd failed: return with ok = .false., so that the caller
! can use a more robust algorithm.
!
            if (it .gt. maxit) go to 100
            call dgemm('t', 'n', m, m, n, one, u, n, u, n, zero, metric, m)
!
! norms of the vectors, and scaling of the metric. vectors that vanish (or contain NaN)
! cannot be orthonormalized: return with an error status.
!
            do j = 1, m
                norms(j) = sqrt(metric(j, j))
            end do
            if (it .eq. 1 .and. present(col_norms)) col_norms = norms
            if (.not. all(norms .gt. zero)) go to 100
            do j = 1, m
                do i = 1, m
                    metric(i, j) = metric(i, j)/(norms(i)*norms(j))
                end do
            end do
            msave = metric
!
! compute the cholesky factorization of the metric.
!
            call dpotrf('l', m, metric, m, info)
!
! if dpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!
            shifted = info .ne. 0
            if (info .ne. 0) then
!
                alpha = 100.0_dp
                it_micro = 0
                micro_done = .false.
!
! add larger and larger shifts to the diagonal until dpotrf manages to factorize it.
!
                do while (.not. micro_done)
                    it_micro = it_micro + 1
!
! something went very wrong. return with an error status, the orthogonalization
! will be carried out using a different algorithm.
!
                    if (it_micro .gt. maxit) go to 100
!
                    shift = max(epsilon(one)*alpha, tol_ortho)
                    metric = msave
                    call diag_shift(m, shift, metric)
                    call dpotrf('l', m, metric, m, info)
                    alpha = alpha*10.0_dp
                    micro_done = info .eq. 0
                end do
!
            end if
!
! we assume that the error on the orthogonality is of order k(l)^2 * eps,
! where eps is the machine precision.
! the condition number k(l) is estimated by computing
!
! k(l) ||l|| ||l^-1||,
!
! where the norm used is the following (see norm_estimate):
!
! || A || = || D + O || <= || D ||_inf + || O ||_2
!
! compute l^-1, using msave to store the inverse cholesky factor
!
            msave = metric
            call dtrtri('l', 'n', m, msave, m, info)
!
! compute the norm of l, l^-1 and the condition number:
!
            l_norm = norm_est(m, metric)
            linv_norm = norm_est(m, msave)
            rcond = l_norm*linv_norm
!
! in each iteration of ortho_cd, we apply l^-t to u, which introduces
! a numerical error of order ||l^-1||.
! this error is saved in growth and used in ortho_vs_x to check how much
! ortho_cd spoiled the previously computed orthogonality to x.
!
            if (present(max_growth)) then
                if (.not. linv_norm .le. max_growth) go to 100
            end if
            growth = growth*linv_norm
!
! orthogonalize u by applying (l^-1 d^-1)^t = d^-1 l^-t, where d holds the norms
!
            do j = 1, m
                msave(:, j) = msave(:, j)/norms(j)
            end do
            call dtrmm('r', 'l', 't', 'n', n, m, one, msave, m, u, n)
!
! check the error:
!
            error = epsilon(one)*rcond*rcond
!
! an iteration that needed a level shift does not orthonormalize the vectors (e.g., zero
! vectors stay zero): another one is always needed.
!
            macro_done = error .lt. tol_ortho .and. .not. shifted
        end do
!
        ok = .true.
!
100     continue
        call mfree(ctx, metric)
        call mfree(ctx, msave)
        call mfree(ctx, norms)
!
    end subroutine ortho_cd
!
    subroutine biortho_eigvecs(ctx, n, m, v_l, v_r)
!* Biorthonormalize a set of left and right eigenvectors, so that
! \[ v_l^Tv_r = \textbf{I}. \]
!
! Left and right eigenvectors that belong to different eigenvalues are already orthogonal
! (up to the convergence error), so they must not be mixed, as, e.g., an SVD of the overlap
! would do. The overlap is instead factorized as
! \[ v_l^Tv_r = LU, \]
! with \(L\) lower triangular and \(U\) unit upper triangular, and
! \[ v_r = v_r U^{-1}, \quad v_l = v_l L^{-T}. \]
! The off-diagonal elements of the factors are of the order of the convergence error, unless
! some eigenvalues are degenerate: in this case, the vectors of the degenerate subspace are
! also biorthogonalized among themselves. Finally, the right eigenvectors are normalized, and
! the left ones scaled accordingly.
!
        implicit none
        type(dgl_context), intent(inout) :: ctx
!
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of vectors
        real(dp), dimension(n, m), intent(inout) :: v_l
!! Left eigenvectors
        real(dp), dimension(n, m), intent(inout) :: v_r
!! Right eigenvectors
!
        integer(ip) :: i, j
        real(dp) :: thresh, fac
        real(dp), allocatable :: over(:, :)
!
        call mallocate(ctx, m, m, over)
        if (dgl_failed(ctx)) go to 100
!
! compute the overlap and factorize it in place (crout, without pivoting, which would
! reorder the eigenvectors):
!
        call dgemm('t', 'n', m, m, n, one, v_l, n, v_r, n, zero, over, m)
!
        do j = 1, m
            do i = j, m
                over(i, j) = over(i, j) - dot_product(over(i, 1:j - 1), over(1:j - 1, j))
            end do
!
! a vanishing pivot means that the left and right eigenvectors are (numerically) orthogonal,
! i.e., that the eigenvalue is defective: they cannot be biorthonormalized.
!
            thresh = 1.0e2_dp*epsilon(one)*dnrm2(n, v_l(1, j), 1_ip)*dnrm2(n, v_r(1, j), 1_ip)
            if (abs(over(j, j)) .le. thresh) then
                call dgl_error(ctx, 'biortho_eigvecs: left and right eigenvectors are orthogonal.', dgl_err_ortho)
                go to 100
            end if
!
            do i = j + 1, m
                over(j, i) = (over(j, i) - dot_product(over(j, 1:j - 1), over(1:j - 1, i)))/over(j, j)
            end do
        end do
!
! v_r = v_r U^-1, v_l = v_l L^-T
!
        call dtrsm('r', 'u', 'n', 'u', n, m, one, over, m, v_r, n)
        call dtrsm('r', 'l', 't', 'n', n, m, one, over, m, v_l, n)
!
        do j = 1, m
            fac = dnrm2(n, v_r(1, j), 1_ip)
            v_r(:, j) = v_r(:, j)/fac
            v_l(:, j) = v_l(:, j)*fac
        end do
!
100     continue
        call mfree(ctx, over)
    end subroutine biortho_eigvecs
!
    real(dp) function norm_est(m, a)
!* Compute a cheap estimate of the norm of a lower triangular matrix.
! Let \(a = d + o\), where \(d = diag(a)\). Since:
! \[
! || a || \leq || d || + || o ||
! \]
! We compute \(|| d ||\) as \(max_i |d(i)|\) and \(|| o ||\) as its frobenius norm.
!
! This is tight enough, and goes to 1 when \(a\) approaches the identity.
!
        implicit none
        integer(ip), intent(in) :: m
!! Dimension of the matrix
        real(dp), dimension(m, m), intent(in) :: a
!! Matrix to compute the norm
!
! Local vars
!
        integer(ip) :: i, j
        real(dp) :: diag_norm, od_norm
!
        diag_norm = zero
        do i = 1, m
            diag_norm = max(diag_norm, abs(a(i, i)))
        end do
!
        od_norm = zero
        do i = 1, m
            do j = 1, i - 1
                od_norm = od_norm + a(i, j)**2
            end do
        end do
        od_norm = sqrt(od_norm)
!
        norm_est = diag_norm + od_norm
        return
    end function norm_est
!
    subroutine gs_fallback(ctx, n, m, k, x, bx, u)
!* Fallback of the orthonormalization routines, used when the Cholesky-based orthonormalization
! fails or would amplify the vectors too much (see max_growth), i.e., when the vectors are
! (numerically) linearly dependent.
!
! Orthonormalize the vectors \(u(n,k)\) one at a time with Gram-Schmidt, projecting each of
! them against \(x(n,m)\), \( u_j = u_j - x (bx^Tu_j) \), and against the previous
! vectors of \(u\). The vectors that are numerically linearly dependent are replaced with random
! vectors, which are orthonormalized in the same way.
! This is slower than the default route for long vectors, as it uses level-2 BLAS.
!
! Linearly dependent vectors, e.g., preconditioned residuals that span fewer directions than
! their number, cannot be orthonormalized as a block: the orthonormalization amplifies the
! round-off noise, which then is no longer orthogonal to \(x\).
! As in the "twice is enough" algorithm (Kahan, Parlett), if the projection removes most of the
! norm of a vector, the projection is repeated, which makes the result orthogonal to the other
! vectors up to round-off errors of the order of machine precision times its norm before the
! second projection. The vector is linearly dependent if the second projection leaves nothing
! but such errors, or if the projections left nothing but round-off of its original norm.
! Small parts of a vector that survive the projections are kept, as they are mostly made of
! the round-off errors of the first projection, which are distributed like the entries of the
! original vector, and are more useful to expand the space than a random vector.
!
        implicit none
        type(dgl_context), intent(inout) :: ctx
!
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of reference vectors
        integer(ip), intent(in) :: k
!! Number of vectors to orthonormalize
        real(dp), dimension(n, m), intent(in) :: x
!! Reference vectors
        real(dp), dimension(n, m), intent(in) :: bx
!! Vectors that define the projector: \(x\), or \(Bx\) for a B-orthogonalization
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthonormalize
!
        integer(ip) :: j, pass, n_random
        real(dp) :: u_norm, prev_norm, start_norm
        logical :: accepted
        real(dp), allocatable :: xu(:), uu(:)
!
        integer(ip), parameter :: maxit = 10
!
        call mallocate(ctx, m, xu)
        call mallocate(ctx, k, uu)
        if (dgl_failed(ctx)) go to 100
!
        do j = 1, k
            n_random = 0
            do
                prev_norm = dnrm2(n, u(1, j), 1_ip)
                if (ieee_is_nan(prev_norm)) then
                    call dgl_error(ctx, 'the vectors to orthogonalize contain NaN.', dgl_err_ortho)
                    go to 100
                end if
!
                start_norm = prev_norm
                accepted = .false.
                if (prev_norm .gt. zero) then
                    do pass = 1, 2
!
! u_j = u_j - x (bx^t u_j) - u_(1:j-1) (u_(1:j-1)^t u_j)
!
                        call dgemv('t', n, m, one, bx, n, u(1, j), 1_ip, zero, xu, 1_ip)
                        call dgemv('n', n, m, -one, x, n, xu, 1_ip, one, u(1, j), 1_ip)
                        if (j .gt. 1) then
                            call dgemv('t', n, j - 1, one, u, n, u(1, j), 1_ip, zero, uu, 1_ip)
                            call dgemv('n', n, j - 1, -one, u, n, uu, 1_ip, one, u(1, j), 1_ip)
                        end if
                        u_norm = dnrm2(n, u(1, j), 1_ip)
                        if (pass .eq. 1) then
                            accepted = u_norm .ge. keep_ratio*prev_norm
                        else
                            accepted = u_norm .gt. real(m + j, dp)*epsilon(one)*prev_norm
                        end if
!
! a vector of which the projections left nothing but round-off of its original norm is
! linearly dependent on x (or on the previous vectors): what is left is round-off with
! components along them, which normalizing would turn into a vector inside their span
!
                        if (u_norm .le. real(m + k, dp)*epsilon(one)*start_norm) accepted = .false.
                        if (accepted .and. u_norm .gt. zero) exit
                        accepted = .false.
                        prev_norm = u_norm
                    end do
                end if
!
                if (accepted) exit
!
! u_j is linearly dependent: replace it with a random vector.
!
                n_random = n_random + 1
                if (n_random .gt. maxit) then
                    call dgl_error(ctx, 'unable to replace linearly dependent vectors.', dgl_err_ortho)
                    go to 100
                end if
                call random_number(u(:, j))
                u(:, j) = u(:, j) - 0.5_dp
            end do
            u(:, j) = u(:, j)/u_norm
        end do
!
100     continue
        call mfree(ctx, xu)
        call mfree(ctx, uu)
    end subroutine gs_fallback
!
    subroutine ortho_gs(ctx, n, k, u)
!* Orthonormalize the vectors \(u(n,k)\), e.g., guess vectors provided by the user, which can be
! incomplete (e.g., zero or repeated vectors). The Cholesky-based orthonormalization is used,
! with the Gram-Schmidt fallback, which replaces the linearly dependent vectors with random
! ones, if it fails.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: k
!! Number of vectors
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthonormalize
!
        real(dp) :: no_x(n, 0)
!
        call orthonormalize_vs_x(ctx, n, 0_ip, k, no_x, no_x, u, "ortho_gs")
    end subroutine ortho_gs
!
    subroutine ortho_vs_x(ctx, n, m, k, x, u)
!* Given two sets \(x(n,m)\) and \(u(n,k)\) of vectors, where \(x\)
! is assumed to be orthonormal, orthogonalize \(u\) against \(x\), and orthonormalize \(u\).
! See orthonormalize_vs_x.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of reference vectors
        integer(ip), intent(in) :: k
!! Number of vectors to orthogonalize against \(x\)
        real(dp), dimension(n, m), intent(in) :: x
!! Reference vectors
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthogonalize
!
        call orthonormalize_vs_x(ctx, n, m, k, x, x, u, "ortho_vs_x")
    end subroutine ortho_vs_x
!
    subroutine b_ortho_vs_x(ctx, n, m, k, x, bx, u)
!*  Given two sets \(x(n,m)\) and \(u(n,k)\) of vectors, where \(x\)
! is assumed to be B-orthonormal, B-orthogonalize \(u\) against \(x\), and orthonormalize \(u\).
! See orthonormalize_vs_x.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of reference vectors
        integer(ip), intent(in) :: k
!! Number of vectors to orthogonalize
        real(dp), dimension(n, m), intent(in) :: x
!! Reference vectors
        real(dp), dimension(n, m), intent(in) :: bx
!! Application of an external \(B\) matrix to \(x\)
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthogonalize
!
        call orthonormalize_vs_x(ctx, n, m, k, x, bx, u, "b_ortho_vs_x")
    end subroutine b_ortho_vs_x
!
    subroutine orthonormalize_vs_x(ctx, n, m, k, x, bx, u, caller)
!* Orthogonalize \(u(n,k)\) against \(x(n,m)\), \( u = u - x (bx^Tu) \), and orthonormalize \(u\),
! iteratively, until the overlap between \(x\) and the orthonormalized \(u\) is smaller than a
! (tight) threshold. With \(m = 0\), just orthonormalize \(u\).
!
! The vectors are first orthonormalized (ortho_cd), which makes them stay nearly orthonormal
! after the projection, so that a single projection is usually enough.
! Each iteration then projects \(u\) against \(x\) (level-3 BLAS) and orthonormalizes it again
! with ortho_cd, which also returns the norms the vectors had after the projection.
! The orthogonality error with respect to \(x\) is estimated as machine precision times the
! growth factor, i.e., the norm of the transformation applied to \(u\) after the projection:
! the inverse of the norms of the projected vectors, times the growth factor of ortho_cd.
!
! If the vectors are (numerically) linearly dependent, the gram-schmidt fallback (gs_fallback)
! is used, which also replaces the dependent vectors with random vectors. This happens if:
! - a vector vanishes after the projection;
! - ortho_cd fails, e.g., because one of its transformations, which only depend on the
!   conditioning of the normalized vectors, has a norm larger than max_growth_cd;
! - the projection leaves nothing but round-off of a vector, i.e., the vector is linearly
!   dependent on \(x\);
! - the default route does not converge within max_default_it iterations.
!
        implicit none
        type(dgl_context), intent(inout) :: ctx
!
        integer(ip), intent(in) :: n
!! Length of the vectors
        integer(ip), intent(in) :: m
!! Number of reference vectors
        integer(ip), intent(in) :: k
!! Number of vectors to orthogonalize
        real(dp), dimension(n, m), intent(in) :: x
!! Reference vectors
        real(dp), dimension(n, m), intent(in) :: bx
!! \(x\), or the application of an external \(B\) matrix to \(x\)
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthogonalize
        character(len=*), intent(in) :: caller
!! Name of the calling routine, for error messages
!
! local variables:
! ================
!
        logical :: done, ok, fallback, used_fallback
        integer(ip) :: it
        real(dp) :: growth, cd_growth
        real(dp), allocatable :: xu(:, :), norm_old(:), norm_new(:)
!
        integer(ip), parameter :: maxit = 10
!
        call mallocate(ctx, m, k, xu)
        call mallocate(ctx, k, norm_old)
        call mallocate(ctx, k, norm_new)
        if (dgl_failed(ctx)) go to 100
!
! orthonormalize the vectors: after the projection, they then stay nearly orthonormal, and
! one projection is usually enough. ortho_cd also returns the norms of the vectors in input.
!
        call ortho_cd(ctx, n, k, u, cd_growth, ok, max_growth_cd, norm_old)
        if (.not. ok) call gs_fallback(ctx, n, m, k, x, bx, u)
        if (dgl_failed(ctx)) go to 100
        norm_old = one
        done = m .eq. 0
        it = 0
        used_fallback = .false.
!
        do while (.not. done)
            it = it + 1
!
! if things went really wrong, abort.
!
            if (it .gt. maxit) then
                call dgl_error(ctx, 'catastrophic failure of '//caller, dgl_err_ortho)
                go to 100
            end if
!
! u = u - x (bx^t u)
!
            call dgemm('t', 'n', m, k, n, one, bx, n, u, n, zero, xu, m)
            call dgemm('n', 'n', n, k, m, -one, x, n, xu, m, one, u, n)
!
! orthonormalize the projected vectors: the norms they had before being orthonormalized,
! returned in norm_new, are needed to estimate the orthogonality error, and to detect the
! vectors that are linearly dependent on x
!
            call ortho_cd(ctx, n, k, u, cd_growth, ok, max_growth_cd, norm_new)
            if (dgl_failed(ctx)) go to 100
            fallback = .not. ok
            if (ok) fallback = any(norm_new .le. real(m + k, dp)*epsilon(one)*norm_old)
!
! the vectors are (numerically) linearly dependent: use the gram-schmidt fallback, which
! returns orthonormal vectors, and check the result in the next iteration.
!
            if (fallback) then
                call gs_fallback(ctx, n, m, k, x, bx, u)
                if (dgl_failed(ctx)) go to 100
                used_fallback = .true.
                cycle
            end if
!
! the orthonormalization has introduced an error that makes the new vectors no longer fully
! orthogonal to x. assuming that u was orthogonal to x to machine precision after the
! projection, the error is estimated as growth * eps, where growth accounts for the scaling of
! the projected vectors and for the transformations applied by ortho_cd.
!
            growth = cd_growth*maxval(norm_old/norm_new)
            if (ieee_is_nan(growth)) then
                call dgl_error(ctx, caller//': the orthogonalization produced NaN.', dgl_err_ortho)
                go to 100
            end if
            done = growth*epsilon(one) .lt. tol_ortho
!
! the default route is not converging: the vectors are (numerically) linearly dependent,
! even if the checks above did not detect it (e.g., with a B-metric, where the projector is
! not orthogonal and the norm of a dependent vector does not vanish).
!
            if (.not. done .and. it .ge. max_default_it .and. .not. used_fallback) then
                call gs_fallback(ctx, n, m, k, x, bx, u)
                if (dgl_failed(ctx)) go to 100
                used_fallback = .true.
            end if
        end do
!
100     continue
        call mfree(ctx, xu)
        call mfree(ctx, norm_old)
        call mfree(ctx, norm_new)
    end subroutine orthonormalize_vs_x

end module dgl_orthogonalizations
