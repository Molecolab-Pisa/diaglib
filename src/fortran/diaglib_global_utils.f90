module diaglib_global_utils
  implicit none
!
! global kind for double precision
!
integer, parameter :: dp = selected_real_kind(15)
!
! useful constants
!
  real(dp), parameter    :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
!
! convergence thresholds for orthogonalization
!
  real(dp), parameter    :: tol_ortho = two * epsilon(one)
! 
! memory and info for lapack routines
!
  integer                :: lwork, info
  real(dp), allocatable  :: work(:), tau(:)
!
! timings:
!
  real(dp)               :: t1(2), t2(2), t_diag(2), t_ortho(2), &
                            t_mv(2), t_tot1(2), t_tot2(2), t_tot(2)
!
!
!
  contains
!
  subroutine check_mem(istat)
    implicit none
!
!   silly routine that checks that istat is zero after memory allocation
!   and aborts otherwise printing an error message.
!
    integer, intent(in) :: istat
!
    1000 format(t3,'memory allocation failed. istat = ',i8)
    if (istat.ne.0) then
      write(6,1000) istat
      stop
    end if
    return
  end subroutine check_mem
!
  integer function get_mem_lapack(n,n_max)
    integer, intent(in)    :: n, n_max
!
    integer           :: lwork1, lwork2, len_rr, len_qr, nb
!fl
    integer           :: lwork3
    integer, external :: ilaenv
!
!   maximum size of the rayleigh-ritz matrix:
!
    len_rr = 3*n_max
!
!   maximum size of the space to orthonormalize:
!
    len_qr = 6*n_max
!
!   use lapack query routines to compute the optimal memory required
!   for diagonalization and QR decomposition.
!
    nb     = ilaenv( 1, 'DSYTRD', 'l', len_rr, -1, -1, -1 )
    lwork1 = len_rr * nb
!
    nb     = ilaenv( 1, 'DGEQRF', 'l', n, len_qr, -1, -1 )
    lwork2 = len_qr*nb
!
    nb     = ilaenv( 1, 'DSYTRD', 'l', len_rr, -1, -1, -1)
    lwork3 = len_rr * nb
!
    get_mem_lapack = max(lwork1,lwork2,lwork3)
    return
  end function get_mem_lapack
! 
  subroutine prtmat(n,lda,mat,iform)
    integer,                      intent(in) :: n, lda, iform
    real(dp), dimension(lda,lda), intent(in) :: mat
!
    integer :: i
!
    100 format(t3,20f12.6)
    200 format(t3,20d12.4)
!
    do i = 1, n
      if (iform.eq.1) write(6,100) mat(i,1:n)
      if (iform.eq.2) write(6,200) mat(i,1:n)
    end do
    return
  end subroutine prtmat
!
  subroutine get_time(t)
    real(dp), dimension(2), intent(inout) :: t
!
!$  real(dp) :: omp_get_wtime
!$  external :: omp_get_wtime
!
!   get cpu and (if openmp is available) wall time.
!
    t = zero
    call cpu_time(t(1))
!$  t(2) = omp_get_wtime()
!
    return
  end subroutine get_time

end module diaglib_global_utils
