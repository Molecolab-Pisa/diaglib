module dgl_lapack
!* Explicit interfaces to the BLAS and LAPACK routines used in DiagLib.
! All the integer arguments are integer(ip), i.e., of the same kind as the integers of the
! BLAS/LAPACK library DiagLib is linked to: with these interfaces, the compiler checks that
! every call passes integers of the right kind (e.g., 1_ip and not 1), so that DiagLib
! does not need to change the default integer kind to be built with 64-bit integers.
    use dgl_interface, only: ip => dgl_int, dp => dgl_real
    implicit none
!
    interface
!
! BLAS
!
        subroutine daxpy(n, da, dx, incx, dy, incy)
            import :: ip, dp
            integer(ip), intent(in) :: n, incx, incy
            real(dp), intent(in) :: da, dx(*)
            real(dp), intent(inout) :: dy(*)
        end subroutine daxpy
!
        subroutine dcopy(n, dx, incx, dy, incy)
            import :: ip, dp
            integer(ip), intent(in) :: n, incx, incy
            real(dp), intent(in) :: dx(*)
            real(dp), intent(inout) :: dy(*)
        end subroutine dcopy
!
        real(dp) function dnrm2(n, x, incx)
            import :: ip, dp
            integer(ip), intent(in) :: n, incx
            real(dp), intent(in) :: x(*)
        end function dnrm2
!
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: ip, dp
            character(len=1), intent(in) :: trans
            integer(ip), intent(in) :: m, n, lda, incx, incy
            real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
            real(dp), intent(inout) :: y(*)
        end subroutine dgemv
!
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: ip, dp
            character(len=1), intent(in) :: transa, transb
            integer(ip), intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
!
        subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            import :: ip, dp
            character(len=1), intent(in) :: side, uplo, transa, diag
            integer(ip), intent(in) :: m, n, lda, ldb
            real(dp), intent(in) :: alpha, a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
        end subroutine dtrsm
!
        subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            import :: ip, dp
            character(len=1), intent(in) :: side, uplo, transa, diag
            integer(ip), intent(in) :: m, n, lda, ldb
            real(dp), intent(in) :: alpha, a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
        end subroutine dtrmm
!
! LAPACK
!
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: ip, dp
            character(len=1), intent(in) :: uplo
            integer(ip), intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer(ip), intent(out) :: info
        end subroutine dpotrf
!
        subroutine dtrtri(uplo, diag, n, a, lda, info)
            import :: ip, dp
            character(len=1), intent(in) :: uplo, diag
            integer(ip), intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer(ip), intent(out) :: info
        end subroutine dtrtri
!
        subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
            import :: ip, dp
            integer(ip), intent(in) :: m, n, lda, lwork
            real(dp), intent(inout) :: a(lda, *), tau(*), work(*)
            integer(ip), intent(out) :: info
        end subroutine dgeqrf
!
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: ip, dp
            character(len=1), intent(in) :: jobz, uplo
            integer(ip), intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *), w(*), work(*)
            integer(ip), intent(out) :: info
        end subroutine dsyev
!
        subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
            import :: ip, dp
            integer(ip), intent(in) :: itype, n, lda, ldb, lwork
            character(len=1), intent(in) :: jobz, uplo
            real(dp), intent(inout) :: a(lda, *), b(ldb, *), w(*), work(*)
            integer(ip), intent(out) :: info
        end subroutine dsygv
!
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            import :: ip, dp
            character(len=1), intent(in) :: jobvl, jobvr
            integer(ip), intent(in) :: n, lda, ldvl, ldvr, lwork
            real(dp), intent(inout) :: a(lda, *), wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer(ip), intent(out) :: info
        end subroutine dgeev
!
        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import :: ip, dp
            character(len=1), intent(in) :: jobu, jobvt
            integer(ip), intent(in) :: m, n, lda, ldu, ldvt, lwork
            real(dp), intent(inout) :: a(lda, *), s(*), u(ldu, *), vt(ldvt, *), work(*)
            integer(ip), intent(out) :: info
        end subroutine dgesvd
!
        integer(ip) function ilaenv(ispec, name, opts, n1, n2, n3, n4)
            import :: ip
            integer(ip), intent(in) :: ispec, n1, n2, n3, n4
            character(len=*), intent(in) :: name, opts
        end function ilaenv
!
    end interface
!
end module dgl_lapack
