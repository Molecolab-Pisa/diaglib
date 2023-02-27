program main
  use real_precision
  use utils
  use diaglib
!
! simple minded test program to call the lobpcg routine
!
  implicit none
  integer  :: n, n_eig, itmax, i, j, n_want, m_max
  real(dp) :: shift, tol
  logical  :: ok
  external :: mmult, mprec
!
  real(dp), allocatable :: eig(:), evec(:,:), diagonal(:)
!
! initialize:
!
  n      = 1000
  n_want = 20
  n_eig  = min(2*n_want,n_want+5)
  tol    = 1.0e-10_dp
  shift  = 0.0_dp
  itmax  = 1000
  nmult  = 0
  m_max  = 25
!
! allocate memory. note that the matrix a is defined in the stack of module "utils".
!
  allocate (a(n,n), eig(n_eig), evec(n,n_eig), diagonal(n))
!
! build a silly matrix:
!
  do i = 1, n
    a(i,i) = dble(i) + 1.0_dp
    diagonal(i) = a(i,i)
    do j = 1, i - 1
      call random_number(a(j,i)) 
      a(i,j) = a(j,i)
    end do
  end do
!
! clean up arrays, make a very crude guess of the first eigenvector:
!
  eig  = 0.0_dp
  evec = 0.0_dp
  do i = 1, n_eig
    evec(i,i) = 1.0_dp
  end do
  call lobpcg_driver(.true.,n,n_want,n_eig,itmax,tol,shift,mmult,mprec,eig,evec,ok)
!
  write(6,*) ' # matmul = ', nmult
  write(6,*)
  nmult = 0
!
! clean up arrays, make a very crude guess of the first eigenvector:
!
  eig  = 0.0_dp
  evec = 0.0_dp
  do i = 1, n_eig
    evec(i,i) = 1.0_dp
  end do
  call davidson_driver(.true.,n,n_want,n_eig,itmax,tol,m_max,shift,mmult,mprec,eig,evec,ok)
!
  write(6,*) ' # matmul = ', nmult
!
end program main
  subroutine mmult(n,m,x,ax)
    use utils
    implicit none
!
!   simple-minded matrix-vector multiplication routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,m), intent(inout) :: ax
!
    integer :: icol
!
    nmult = nmult + m
    do icol = 1, m
      ax(:,icol) = matmul(a,x(:,icol))
    end do
    return
  end subroutine mmult
!
  subroutine mprec(n,m,fac,x,px)
    use utils
    implicit none
!
!   simple-minded shift-and-invert preconditioner routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                   intent(in)    :: n, m
    real(dp),                  intent(in)    :: fac
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,m), intent(inout) :: px
!
    integer             :: i, icol
    real(dp), parameter :: tol = 1.0d-5
!
    do icol = 1, m
      do i = 1, n
        if (abs(a(i,i)+fac).gt.tol) px(i,icol) = x(i,icol)/(a(i,i) + fac)
      end do
    end do
    return
  end subroutine mprec
!
