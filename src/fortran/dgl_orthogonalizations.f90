module dgl_orthogonalizations
!* Module for all the orthogonalization procedures.
! They may be used also as a standalone outside of DiagLib.
    use dgl_global_utils
    implicit none
!
    real(dp), parameter, private :: tol_ortho = two*epsilon(one)
!! Convergence thresholds for orthogonalizations
!
contains
!
    subroutine ortho(n, m, u, w)
!! Orthogonalization routine based on QR decomposition.
!! Orthogonalizes \(m\) vectors of lenght \(n\) contained in \(u\).
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
!
        integer, intent(in) :: n
!! Lenght of the input vectors
        integer, intent(in) :: m
!! Number of input vectors
        real(dp), dimension(n, m), intent(inout) :: u
!! Vectors to orthogonalize
        real(dp), dimension(n, m), optional, intent(inout) :: w
!! Second set of vectors to which we may apply the same transformation
!
!   local scratch
!   =============
!
        real(dp), allocatable :: v(:, :)
!
        call mallocate(n, m, v)
        v = u
        call dgeqrf(n, m, u, n, tau, work, lwork, info)
!
        call dtrsm('r', 'u', 'n', 'n', n, m, one, u, n, v, n)
!
        if (present(w)) call dtrsm('r', 'u', 'n', 'n', n, m, one, u, n, w, n)
!
        u = v
!
        call mfree(v)
    end subroutine ortho
!
    subroutine b_ortho(n, m, u, bu)
!! Subroutine to B-orthogonalize \(m\) vectors of lenght \(n\)
!! using the Cholesky factorization of their overlap.
!! \[ u^TBu = \textbf{I} \]
!!
!! This is in principle not a good idea, as the \(u^TBu\) matrix can be very
!! ill-conditioned, independently of how bad is \(B\). Also, only works if u is
!! already orthonormal.
!! More details on the orthoganlization scheme are available in the documentation
!! of the [[ortho_cd]] procedure
        implicit none
!
        integer, intent(in) :: n
!! Lenght of the vectors
        integer, intent(in) :: m
!! Number of vectors
        real(dp), dimension(n, m), intent(inout) :: u
!! Vectors to B-orthogonalize
        real(dp), dimension(n, m), intent(inout) :: bu
!! Application of an external matrix \(B\) on \(u\)
!
!   local variables
!   ===============
!
        integer :: info, i, j
        real(dp), allocatable :: metric(:, :), sigma(:), u_svd(:, :), vt_svd(:, :), &
                                 temp(:, :)
        real(dp), parameter :: tol_svd = 1.0e-5_dp
        logical, parameter :: use_svd = .false.
!
!
        call mallocate(m, m, metric)
!
        call dgemm('t', 'n', m, m, n, one, u, n, bu, n, zero, metric, m)
!
        if (use_svd) then
!
!     debug option: use svd to b-orthonormalize, by computing
!     b**(-1/2)
!
            call mallocate(m, sigma)
            call mallocate(m, m, u_svd)
            call mallocate(m, m, vt_svd)
            call mallocate(n, m, temp)
!
            call dgesvd('a', 'a', m, m, metric, m, sigma, u_svd, m, vt_svd, m, work, lwork, info)
!
!     compute sigma**(-1/2)
!
            do i = 1, m
                if (sigma(i) .gt. tol_svd) then
                    sigma(i) = 1/sqrt(sigma(i))
                else
                    sigma(i) = zero
                end if
            end do
!
!     compute metric ** (-1/2). first, compute sigma ** (-1/2) vt
!
            metric = zero
            do i = 1, m
                do j = 1, m
                    metric(j, i) = metric(j, i) + sigma(j)*vt_svd(j, i)
                end do
            end do
!
!     now, multiply for u:
!
            vt_svd = metric
            call dgemm('n', 'n', m, m, m, one, u_svd, m, vt_svd, m, zero, metric, m)
!
!     metric contains s ** (-1/2), and projects out directions corresponding
!     to pathological singular values.
!     orthogonalize u and bu:
!
            call dgemm('n', 'n', n, m, m, one, u, n, metric, m, zero, temp, n)
            u = temp
            call dgemm('n', 'n', n, m, m, one, bu, n, metric, m, zero, temp, n)
            bu = temp
!
            call mfree(sigma)
            call mfree(u_svd)
            call mfree(vt_svd)
            call mfree(temp)
        else
!
!     compute the cholesky factorization of the metric.
!
            call dpotrf('l', m, metric, m, info)
!
!     get u * l^-T and bu * l^-T
!
            call dtrsm('r', 'l', 't', 'n', n, m, one, metric, m, u, n)
            call dtrsm('r', 'l', 't', 'n', n, m, one, metric, m, bu, n)
        end if
!
        call mfree(metric)
!
    end subroutine b_ortho
!
    subroutine diag_shift(n, shift, a)
!! Add a shift to the diagonal elements of the matrix \(a\)
        implicit none
!
        integer, intent(in) :: n
        real(dp), intent(in) :: shift
        real(dp), dimension(n, n), intent(inout) :: a
!
        integer :: i
!
        do i = 1, n
            a(i, i) = a(i, i) + shift
        end do
!
        return
    end subroutine diag_shift
!
    subroutine ortho_cd(n, m, u, growth, ok)
!! Subroutine to orthogonalize \(m\) vectors of lenght \(n\)
!! using the Cholesky factorization of their overlap.
!! \[ u^Tu = \textbf{I} \]
!! The metric is computed as \(metric = u^Tu \) and then by computing its cholesky
!! decompositoin \( metric = LL^T \). The orthogonal vectors are obtained then
!! by solving the triangular linear system \( u(ortho)L^T = u \).
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
!! more robust routines (QR or SVD) in case of failure.
!
        implicit none
!
        integer, intent(in) :: n
!! Lenght of the vectors
        integer, intent(in) :: m
!! Number of vectors
        real(dp), dimension(n, m), intent(inout) :: u
!! Vectors to orthogonalize
        real(dp), intent(inout) :: growth
!! Growth factor for the numerical error
        logical, intent(inout) :: ok
!! Status of the procedure in output
!
!   local variables
!   ===============
!
        integer :: it, it_micro
        real(dp) :: error, alpha, unorm, shift
        real(dp) :: rcond, l_norm, linv_norm
        logical :: macro_done, micro_done
        real(dp), parameter :: tol_ortho = two*epsilon(one)
        integer, parameter :: maxit = 10
!
!   local scratch
!   =============
!
        real(dp), allocatable :: metric(:, :), msave(:, :)
!
!
!   get memory for the metric.
!
        call mallocate(m, m, metric)
        call mallocate(m, m, msave)
!
        metric = zero
        macro_done = .false.
!
!   assemble the metric
!
        it = 0
        growth = one
        do while (.not. macro_done)
            it = it + 1
            if (it .gt. maxit) then
!
!       ortho_cd failed. return with an error message
!
                ok = .false.
                write (6, 100) ' maximum number of iterations reached.'
                return
            end if
            call dgemm('t', 'n', m, m, n, one, u, n, u, n, zero, metric, m)
            msave = metric
!
!   compute the cholesky factorization of the metric.
!
            call dpotrf('l', m, metric, m, info)
!
!     if dpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!
            if (info .ne. 0) then
!
                alpha = 100.0_dp
                unorm = dnrm2(n*m, u, 1)
                it_micro = 0
                micro_done = .false.
!
!       add larger and larger shifts to the diagonal until dpotrf manages to factorize it.
!
                do while (.not. micro_done)
                    it_micro = it_micro + 1
                    if (it_micro .gt. maxit) then
!
!           something went very wrong. return with an error status, the orthogonalization
!           will be carried out using a different algorithm.
!
                        ok = .false.
                        write (6, 100) ' maximum number of iterations for factorization reached.'
                        stop
                        return
                    end if
!
                    shift = max(epsilon(one)*alpha*unorm, tol_ortho)
                    metric = msave
                    call diag_shift(m, shift, metric)
                    call dpotrf('l', m, metric, m, info)
                    alpha = alpha*10.0_dp
                    micro_done = info .eq. 0
                end do
!
            end if
!
!     we assume that the error on the orthogonality is of order k(l)^2 * eps,
!     where eps is the machine precision.
!     the condition number k(l) is estimated by computing
!
!     k(l) ||l|| ||l^-1||,
!
!     where the norm used is the following (see norm_estimate):
!
!     || A || = || D + O || <= || D ||_inf + || O ||_2
!
!     compute l^-1, using msave to store the inverse cholesky factor
!
            msave = metric
            call dtrtri('l', 'n', m, msave, m, info)
!
!     compute the norm of l, l^-1 and the condition number:
!
            l_norm = norm_est(m, metric)
            linv_norm = norm_est(m, msave)
            rcond = l_norm*linv_norm
!
!     in each iteration of ortho_cd, we apply l^-t to u, which introduces
!     a numerical error of order ||l^-1||.
!     this error is saved in growth and used in ortho_vs_x to check how much
!     ortho_cd spoiled the previously computed orthogonality to x.
!
            growth = growth*linv_norm
!
!     orthogonalize u by applying l^(-t)
!
            call dtrmm('r', 'l', 't', 'n', n, m, one, msave, m, u, n)
!
!     check the error:
!
            error = epsilon(one)*rcond*rcond
            macro_done = error .lt. tol_ortho
        end do
!
100     format(t3, 'ortho_cd failed with the following error:', a)
!
        ok = .true.
!
        call mfree(metric)
        call mfree(msave)
!
    end subroutine ortho_cd
!
    subroutine biortho_vs_x(n, m, k, xl, xr, ul, ur)
!* Given four sets: \(x_l(n,m)\), \(x_r(n,m)\) and \(u_l(n,k)\), \(u_r(n,k)\)
!  of vectors, where \(x_l\) and \(x_r\) are assumed to be orthogonal,
!  orthogonalize \(u_l\) against \(x_l\) and
!  \(u_r\) against \(x_r\).
!
!  Furthermore, orthonormalize \(u_l\) and \(u_r\).
!
        implicit none
!
        integer, intent(in) :: n, m, k
        real(dp), dimension(n, m), intent(in) :: xl, xr
        real(dp), dimension(n, k), intent(inout) :: ul, ur
!
!   local variables:
!
        integer :: it
        real(dp) :: xu_norm(2), growth
        logical :: done, ok
        real(dp), allocatable :: xu(:, :)
!
        integer, parameter :: maxit = 20
!
        call mallocate(m, k, xu)
!
        done = .false.
        it = 0
!
        do while (.not. done)
            it = it + 1
            if (it .gt. maxit) stop 'biortho_vs_x failed.'
!
!     biorthogonalize ul and ur to xr and xl:
!
            call dgemm('t', 'n', m, k, n, one, xl, n, ur, n, zero, xu, m)
            call dgemm('n', 'n', n, k, m, -one, xr, n, xu, m, one, ur, n)
            call dgemm('t', 'n', m, k, n, one, xr, n, ul, n, zero, xu, m)
            call dgemm('n', 'n', n, k, m, -one, xl, n, xu, m, one, ul, n)
!
!     now, orthogonalize ur and ul.
!
            call ortho_cd(n, k, ul, growth, ok)
            xu_norm(1) = growth*epsilon(one)
            call ortho_cd(n, k, ur, growth, ok)
            xu_norm(2) = growth*epsilon(one)
!
            done = xu_norm(1) .lt. tol_ortho .and. xu_norm(2) .lt. tol_ortho
        end do
!
!   make the left and right eigenvectors biorthogonal using the singular value
!   decomposition
!
        call svd_biortho(n, k, ul, ur)
!
        call mfree(xu)
    end subroutine biortho_vs_x
!
    subroutine svd_biortho(n, m, u_l, u_r)
!*  Given two set of vectors, biorthogonalize them by computing the SVD decomposition
!   of the overlap matrix and then solving
! \[
!     metric = u_l^Tu_r \\
!     metric = U \Sigma V^T \\
!     u_l = u_l V^T \\
!     u_r = u_r U
! \]
! Resulting vectors obey:
! \[ u_l^Tu_r = \textbf{I} \]
!
        implicit none
!
        integer, intent(in) :: n
!! Lenght of the vectors
        integer, intent(in) :: m
!! Number of vectors
        real(dp), dimension(n, m), intent(inout) :: u_l
!! First set of vectors
        real(dp), dimension(n, m), intent(inout) :: u_r
!! Second set of vectors
!
        integer :: i
        real(dp) :: fac
!
        real(dp), allocatable :: over(:, :), u(:, :), s(:), vt(:, :), tmp(:, :)
!
!   allocate memory.
!
        call mallocate(m, m, over)
        call mallocate(m, s)
        call mallocate(m, m, u)
        call mallocate(m, m, vt)
        call mallocate(n, m, tmp)
!
!   compute the overlap:
!
        call dgemm('t', 'n', m, m, n, one, u_l, n, u_r, n, zero, over, m)
!
!   compute its singular value decomposition:
!
        call dgesvd('a', 'a', m, m, over, m, s, u, m, vt, m, work, lwork, info)
!
!   compute l*u and r*v
!
        call dgemm('n', 'n', n, m, m, one, u_l, n, u, m, zero, tmp, n)
        u_l = tmp
        call dgemm('n', 't', n, m, m, one, u_r, n, vt, m, zero, tmp, n)
        u_r = tmp
!
!   scale with square root of singular values
!   here, dropping redundant vectors could be a good idea...
!
        do i = 1, m
            fac = one/sqrt(s(i))
            u_l(:, i) = fac*u_l(:, i)
            u_r(:, i) = fac*u_r(:, i)
        end do
!
        call mfree(over)
        call mfree(u)
        call mfree(s)
        call mfree(vt)
        call mfree(tmp)
    end subroutine svd_biortho
!
    real(dp) function norm_est(m, a)
!* Compute a cheap estimate of the norm of a lower triangular matrix.
!  Let \(a = d + o\), where \(d = diag(a)\). Since:
!  \[
!   || a || \leq || d || + || o ||
!  \]
!  We compute \(|| d ||\) as \(max_i |d(i)|\) and \(|| o ||\) as its frobenius norm.
!
!  This is tight enough, and goes to 1 when \(a\) approaches the identity.
!
        implicit none
        integer, intent(in) :: m
!! Dimension of the matrix
        real(dp), dimension(m, m), intent(in) :: a
!! Matrix to compute the norm
!
!   Local vars
!
        integer :: i, j
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
    subroutine ortho_vs_x(n, m, k, x, u, ax, au)
        implicit none
!*  Given two sets \(x(n,m)\) and \(u(n,k)\) of vectors, where \(x\)
!   is assumed to be orthogonal, orthogonalize \(u\) against \(x\).
!
!   If required, orthogonalize au to ax using the same linear
!   transformation, where ax and au are the results of the
!   application of a matrix \(a\) to both \(x\) and \(u\).
!
!   Furthermore, orthonormalize \(u\) and, if required, apply the
!   same transformation to \(au\).
!
!   This routine performs the \(u\) vs \(x\) orthogonalization and the
!   subsequent orthonormalization of \(u\) iteratively, until the
!   overlap between \(x\) and the orthogonalized \(u\) is smaller than
!   a (tight) threshold.
!
        integer, intent(in) :: n
!! Lenght of the vectors
        integer, intent(in) :: m
!! Number of reference vectors
        integer, intent(in) :: k
!! Number of vectors to orthogonalize again \(x\)
        real(dp), dimension(n, m), intent(in) :: x
!! Reference vectors
        real(dp), dimension(n, m), intent(in) :: ax
!! Application of an external matrix on \(x\)
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthogonalize
        real(dp), dimension(n, k), intent(inout) :: au
!! Application of an external matrix on \(u\)
!
!   local variables:
!   ================
!
        logical :: done, ok
        integer :: it
        real(dp) :: xu_norm, growth
        real(dp), allocatable :: xu(:, :)
!
        integer, parameter :: maxit = 10
        logical, parameter :: useqr = .false.
!
!   allocate space for the overlap between x and u.
!
        ok = .false.
        call mallocate(m, k, xu)
        done = .false.
        it = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
        if (.not. useqr) call ortho_cd(n, k, u, growth, ok)
        if (.not. ok .or. useqr) call ortho(n, k, u, au)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
        do while (.not. done)
            it = it + 1
!
!     u = u - x (x^t u)
!
            call dgemm('t', 'n', m, k, n, one, x, n, u, n, zero, xu, m)
            call dgemm('n', 'n', n, k, m, -one, x, n, xu, m, one, u, n)
!
!     now, orthonormalize u.
!
            if (.not. useqr) call ortho_cd(n, k, u, growth, ok)
            if (.not. ok .or. useqr) call ortho(n, k, u, au)
!
!     the orthogonalization has introduced an error that makes the new
!     vector no longer fully orthogonal to x. assuming that u was
!     orthogonal to x to machine precision before, we estimate the
!     error with growth * eps, where growth is the product of the norms
!     of all the linear transformations applied to u.
!     if ortho_cd has failed, we just compute the overlap and its norm.
!
            if (.not. ok .or. useqr) then
                call dgemm('t', 'n', m, k, n, one, x, n, u, n, zero, xu, m)
                xu_norm = dnrm2(m*k, xu, 1)
            else
                xu_norm = growth*epsilon(one)
            end if
            done = xu_norm .lt. tol_ortho
!
!     if things went really wrong, abort.
!
            if (it .gt. maxit) stop ' catastrophic failure of ortho_vs_x'
        end do
!
        call mfree(xu)
!
        return
    end subroutine ortho_vs_x
!
    subroutine b_ortho_vs_x(n, m, k, x, bx, u)
!*  Given two sets \(x(n,m)\) and \(u(n,k)\) of vectors, where \(x\)
!   is assumed to be orthogonal, B-orthogonalize \(u\) against \(x\).
!   furthermore, orthonormalize \(u\).
!
!   This routine performs the \(u\) vs \(x\) orthogonalization and the
!   subsequent orthonormalization of \(u\) iteratively, until the
!   overlap between \(x\) and the orthogonalized \(u\) is smaller than
!   a (tight) threshold.
!
        implicit none
!
        integer, intent(in) :: n
!! Lenght of the vectors
        integer, intent(in) :: m
!! Number of reference vectors
        integer, intent(in) :: k
!! Number of vectors to orthogonalize
        real(dp), dimension(n, m), intent(in) :: x
!! Reference vectors
        real(dp), dimension(n, m), intent(in) :: bx
!! Application of an external \(B\) matrix to \(x\)
        real(dp), dimension(n, k), intent(inout) :: u
!! Vectors to orthogonalize
!
!   local variables:
!   ================
!
        logical :: done, ok
        integer :: it
        real(dp) :: xu_norm, growth, xx(1)
        real(dp), allocatable :: xu(:, :)
!
        integer, parameter :: maxit = 10
        logical, parameter :: useqr = .false.
!
!   allocate space for the overlap between x and u.
!
        ok = .false.
        call mallocate(m, k, xu)
        done = .false.
        it = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
        if (.not. useqr) call ortho_cd(n, k, u, growth, ok)
        if (.not. ok .or. useqr) call ortho(n, k, u)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
        do while (.not. done)
            it = it + 1
!
!     u = u - x (bx^t u)
!
            call dgemm('t', 'n', m, k, n, one, bx, n, u, n, zero, xu, m)
            call dgemm('n', 'n', n, k, m, -one, x, n, xu, m, one, u, n)
!
!     now, orthonormalize u.
!
            if (.not. useqr) call ortho_cd(n, k, u, growth, ok)
            if (.not. ok .or. useqr) call ortho(n, k, u)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
!     note that, if we use ortho_cd, we estimate the norm of the overlap
!     using the growth factor returned in growth.
!     see ortho_vs_x for more information.
!
            if (.not. ok .or. useqr) then
                call dgemm('t', 'n', m, k, n, one, bx, n, u, n, zero, xu, m)
                xu_norm = dnrm2(m*k, xu, 1)
            else
                xu_norm = growth*epsilon(one)
            end if
            done = xu_norm .lt. tol_ortho
!
!     if things went really wrong, abort.
!
            if (it .gt. maxit) stop ' catastrophic failure of b_ortho_vs_x'
        end do
!
        call mfree(xu)
!
    end subroutine b_ortho_vs_x

end module dgl_orthogonalizations
