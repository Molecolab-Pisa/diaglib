module direct_matvecs
    use dgl_interface, only: dgl_real, dgl_int, dgl_success, dgl_err_input, dgl_lobpcg_driver
    use utility, only: close_matrix_element, cplx_matrix_element
!
!   weak_bias is added to the diagonal of dx_close_weak, which then gives a poor
!   preconditioner, so that many iterations and restarts are needed
!
    real(dgl_real), parameter :: weak_bias = 1000.0_dgl_real
!
!   ax_nested runs DiagLib on a smaller problem before applying the matrix: the number
!   of nested calls and of the ones that did not behave as expected are counted here
!
    integer(dgl_int) :: nested_calls = 0, nested_failed = 0
!
contains
!
    subroutine ax(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do i = 1, n
                do j = 1, n
                    if (j .eq. i) then
                        y(i, k) = y(i, k) + real(i + 1, dgl_real)*x(j, k)
                    else
                        y(i, k) = y(i, k) + x(j, k)/real(i + j, dgl_real)
                    end if
                end do
            end do
        end do
!
        return
    end subroutine ax
!
    subroutine mx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do i = 1, n
                do j = 1, n
                    if (j .eq. i) then
                        y(i, k) = y(i, k) + x(j, k)
                    else
                        y(i, k) = y(i, k) + x(j, k)/real(i + j, dgl_real)
                    end if
                end do
            end do
        end do
!
        return
    end subroutine mx
!
    subroutine dx(n, m, shift, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), intent(in) :: shift
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, k
        real(dgl_real) :: fac
!
        real(dgl_real), parameter :: eps = 1.0e-5_dgl_real
!
        do k = 1, m
            do i = 1, n
                fac = shift + real(i + 1, dgl_real)
                if (abs(fac) .gt. eps) then
                    y(i, k) = x(i, k)/(shift + real(i + 1, dgl_real))
                else
                    y(i, k) = x(i, k)
                end if
            end do
        end do
    end subroutine dx
!
    subroutine arx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
        real(dgl_real) :: fac
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do i = 1, n
                do j = 1, n
                    if (j .eq. i) then
                        y(i, k) = y(i, k) + real(i + 1, dgl_real)*x(j, k)
                    else
                        fac = real(i, dgl_real)/real(j, dgl_real)
                        y(i, k) = y(i, k) + fac*x(j, k)/real(i + j, dgl_real)
                    end if
                end do
            end do
        end do
!
        return
    end subroutine arx
!
    subroutine alx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
        real(dgl_real) :: fac
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do i = 1, n
                do j = 1, n
                    if (j .eq. i) then
                        y(i, k) = y(i, k) + real(i + 1, dgl_real)*x(j, k)
                    else
                        fac = real(j, dgl_real)/real(i, dgl_real)
                        y(i, k) = y(i, k) + fac*x(j, k)/real(i + j, dgl_real)
                    end if
                end do
            end do
        end do
!
        return
    end subroutine alx
!
    subroutine arx_close(n, m, x, y)
!
!   non-symmetric matrix with nearly degenerate eigenvalues (see close_matrix_element)
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
        y = 0.0_dgl_real
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    y(i, k) = y(i, k) + close_matrix_element(i, j)*x(j, k)
                end do
            end do
        end do
    end subroutine arx_close
!
    subroutine alx_close(n, m, x, y)
!
!   transpose of the matrix applied by arx_close
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
        y = 0.0_dgl_real
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    y(i, k) = y(i, k) + close_matrix_element(j, i)*x(j, k)
                end do
            end do
        end do
    end subroutine alx_close
!
    subroutine dx_close(n, m, shift, x, y)
!
!   diagonal preconditioner for the matrix applied by arx_close
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), intent(in) :: shift
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, k
        real(dgl_real) :: fac
        real(dgl_real), parameter :: eps = 1.0e-5_dgl_real
!
        do k = 1, m
            do i = 1, n
                fac = shift + close_matrix_element(i, i)
                if (abs(fac) .gt. eps) then
                    y(i, k) = x(i, k)/fac
                else
                    y(i, k) = x(i, k)
                end if
            end do
        end do
    end subroutine dx_close
!
    subroutine dx_close_weak(n, m, shift, x, y)
!
!   poor diagonal preconditioner for the matrix applied by arx_close
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), intent(in) :: shift
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, k
!
        do k = 1, m
            do i = 1, n
                y(i, k) = x(i, k)/(shift + close_matrix_element(i, i) + weak_bias)
            end do
        end do
    end subroutine dx_close_weak
!
    subroutine arx_cplx(n, m, x, y)
!
!   non-symmetric matrix with a pair of complex eigenvalues (see cplx_matrix_element)
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
        y = 0.0_dgl_real
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    y(i, k) = y(i, k) + cplx_matrix_element(i, j)*x(j, k)
                end do
            end do
        end do
    end subroutine arx_cplx
!
    subroutine alx_cplx(n, m, x, y)
!
!   transpose of the matrix applied by arx_cplx
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
        y = 0.0_dgl_real
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    y(i, k) = y(i, k) + cplx_matrix_element(j, i)*x(j, k)
                end do
            end do
        end do
    end subroutine alx_cplx
!
    subroutine sx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
!   just the identity matrix.
!
        y = x
!
        return
    end subroutine sx
!
    subroutine apbx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
!   (a + b)_ij = (5 + i) \delta_ij + (1 - \delta_ij) / (i+j)
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    if (i .eq. j) then
                        y(i, k) = y(i, k) + real(5 + i, dgl_real)*x(i, k)
                    else
                        y(i, k) = y(i, k) + x(j, k)/real(i + j, dgl_real)
                    end if
                end do
            end do
        end do
        return
    end subroutine apbx
!
    subroutine ambx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
!   (a + b)_ij = (2 + i) \delta_ij + (0.2 - \delta_ij) / (i+j)
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    if (i .eq. j) then
                        y(i, k) = y(i, k) + real(2 + i, dgl_real)*x(i, k)
                    else
                        y(i, k) = y(i, k) + 0.20_dgl_real*x(j, k)/real(i + j, dgl_real)
                    end if
                end do
            end do
        end do
        return
    end subroutine ambx
!
    subroutine spdx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
!   \sigma = 1, \delta_ij = +- 0.05
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    if (i .eq. j) then
                        y(i, k) = y(i, k) + x(i, k)
                    else if (i .gt. j) then
                        y(i, k) = y(i, k) + 0.05_dgl_real*x(j, k)
                    else
                        y(i, k) = y(i, k) - 0.05_dgl_real*x(j, k)
                    end if
                end do
            end do
        end do
        return
    end subroutine spdx
!
    subroutine smdx(n, m, x, y)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int) :: i, j, k
!
!   (a + b)_ij = (2 + i) \delta_ij + (0.2 - \delta_ij) / (i+j)
!
        y = 0.0_dgl_real
!
        do k = 1, m
            do j = 1, n
                do i = 1, n
                    if (i .eq. j) then
                        y(i, k) = y(i, k) + x(i, k)
                    else if (i .gt. j) then
                        y(i, k) = y(i, k) - 0.05_dgl_real*x(j, k)
                    else
                        y(i, k) = y(i, k) + 0.05_dgl_real*x(j, k)
                    end if
                end do
            end do
        end do
        return
    end subroutine smdx
!
    subroutine lrprc(n, m, fac, xp, xm, yp, ym)
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), intent(in) :: fac
        real(dgl_real), dimension(n, m), intent(in) :: xp, xm
        real(dgl_real), dimension(n, m), intent(inout) :: yp, ym
!
        integer(dgl_int) :: i, k
        real(dgl_real) :: val
!
!   yp = xp
!   ym = xm
        do k = 1, m
            do i = 1, n
                val = fac*fac*(real(i + 7, dgl_real)**2 - 1.0_dgl_real)
                val = 1.0_dgl_real/val
                yp(i, k) = val*(fac*real(i + 7, dgl_real)*xp(i, k) + xm(i, k))
                ym(i, k) = val*(fac*real(i + 7, dgl_real)*xm(i, k) + xp(i, k))
            end do
        end do
!
        return
    end subroutine lrprc

    real(dgl_real) function rot_eigenvalue(i)
!
!   eigenvalues of ax_rot, in increasing order
!
        implicit none
        integer(dgl_int), intent(in) :: i
        rot_eigenvalue = real(i, dgl_real) + 0.5_dgl_real*sin(real(i, dgl_real))**2
    end function rot_eigenvalue
!
    subroutine householder(n, m, x)
!
!   apply the Householder reflection H = I - 2 v v^T, v_i ~ sin(i), to m vectors
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(inout) :: x
        real(dgl_real) :: v(n)
        integer(dgl_int) :: i, k
!
        do i = 1, n
            v(i) = sin(real(i, dgl_real))
        end do
        v = v/norm2(v)
        do k = 1, m
            x(:, k) = x(:, k) - 2.0_dgl_real*dot_product(v, x(:, k))*v
        end do
    end subroutine householder
!
    subroutine ax_rot(n, m, x, y)
!
!   a = H diag(d) H, with d_i = rot_eigenvalue(i): the eigenvectors are H e_i
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
        integer(dgl_int) :: i
!
        y = x
        call householder(n, m, y)
        do i = 1, n
            y(i, :) = rot_eigenvalue(i)*y(i, :)
        end do
        call householder(n, m, y)
    end subroutine ax_rot
!
    subroutine dx_rot(n, m, shift, x, y)
!
!   diagonal preconditioner for ax_rot, (d + shift)^-1, approximate as the diagonal
!   of ax_rot is not d
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), intent(in) :: shift
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
        integer(dgl_int) :: i
        real(dgl_real) :: den
!
        do i = 1, n
            den = rot_eigenvalue(i) + shift
            if (abs(den) .lt. 1.0e-3_dgl_real) den = sign(1.0e-3_dgl_real, den)
            y(i, :) = x(i, :)/den
        end do
    end subroutine dx_rot
!
    subroutine ax_nan(n, m, x, y)
!
!   a broken matrix-vector product, which returns a NaN
!
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        y = x
        y(n, m) = ieee_value(1.0_dgl_real, ieee_quiet_nan)
    end subroutine ax_nan
!
    subroutine ax_tri(n, m, x, y)
!
!   tridiagonal matrix a_ii = i, a_i,i+1 = a_i+1,i = 0.1. with unit vectors as a guess, the
!   residuals are all proportional to the same unit vector: the preconditioned residuals are
!   linearly dependent
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
        integer(dgl_int) :: i
!
        do i = 1, n
            y(i, :) = real(i, dgl_real)*x(i, :)
        end do
        y(1:n - 1, :) = y(1:n - 1, :) + 0.1_dgl_real*x(2:n, :)
        y(2:n, :) = y(2:n, :) + 0.1_dgl_real*x(1:n - 1, :)
    end subroutine ax_tri
!
    subroutine mx_tri(n, m, x, y)
!
!   tridiagonal metric b_ii = 1, b_i,i+1 = b_i+1,i = 0.05
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        y = x
        y(1:n - 1, :) = y(1:n - 1, :) + 0.05_dgl_real*x(2:n, :)
        y(2:n, :) = y(2:n, :) + 0.05_dgl_real*x(1:n - 1, :)
    end subroutine mx_tri
!
    subroutine dx_tri(n, m, shift, x, y)
!
!   diagonal preconditioner for ax_tri
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), intent(in) :: shift
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
        integer(dgl_int) :: i
        real(dgl_real) :: den
!
        do i = 1, n
            den = real(i, dgl_real) + shift
            if (abs(den) .lt. 1.0e-8_dgl_real) den = 1.0e-8_dgl_real
            y(i, :) = x(i, :)/den
        end do
    end subroutine dx_tri
!
    subroutine ax_nested(n, m, x, y)
!
!   same as ax, but DiagLib is called again from inside the callback, on a smaller
!   linear response matrix: first with invalid input, then to actually solve the problem.
!   neither the error nor the successful inner call may affect the outer call.
!
        implicit none
        integer(dgl_int), intent(in) :: n, m
        real(dgl_real), dimension(n, m), intent(in) :: x
        real(dgl_real), dimension(n, m), intent(inout) :: y
!
        integer(dgl_int), parameter :: n_in = 100, n_targ_in = 2, n_max_in = 4
        real(dgl_real) :: eig(n_max_in), evec(n_in, n_max_in), res(n_in)
        integer(dgl_int) :: i, info
        logical :: ok, passed
!
        nested_calls = nested_calls + 1
        evec = 0.0_dgl_real
        do i = 1, n_max_in
            evec(i, i) = 1.0_dgl_real
        end do
!
        call dgl_lobpcg_driver(n_in, n_max_in + 1, n_max_in, apbx, dx, eig, evec, ok, dgl_info=info)
        passed = info .eq. dgl_err_input .and. .not. ok
!
        call dgl_lobpcg_driver(n_in, n_targ_in, n_max_in, apbx, dx, eig, evec, ok, &
                               dgl_tol=1.0e-8_dgl_real, dgl_info=info)
        passed = passed .and. info .eq. dgl_success .and. ok
        if (passed) then
!
!   check the residual of the first eigenvector
!
            call apbx(n_in, 1_dgl_int, evec(:, 1), res)
            res = res - eig(1)*evec(:, 1)
            passed = norm2(res) .lt. 1.0e-6_dgl_real
        end if
        if (.not. passed) nested_failed = nested_failed + 1
!
        call ax(n, m, x, y)
    end subroutine ax_nested
!
end module direct_matvecs
