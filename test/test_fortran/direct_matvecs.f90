module direct_matvecs
    use dgl_interface, only: dgl_real, dgl_int
    use utility, only: close_matrix_element, cplx_matrix_element
!
!   weak_bias is added to the diagonal of dx_close_weak, which then gives a poor
!   preconditioner, so that many iterations and restarts are needed
!
    real(dgl_real), parameter :: weak_bias = 1000.0_dgl_real
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

end module direct_matvecs
