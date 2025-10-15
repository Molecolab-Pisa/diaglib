program test_fortran
  use real_precision
  use diaglib
  implicit none
!
! tests the various functionalities of diaglib.
!
  integer,  parameter :: n = 1000, n_targ = 20, n_max = 25, max_iter = 100, max_dav = 20
  real(dp), parameter :: tol = 1.0e-8_dp, shift = 0.0_dp, zero = 0.0_dp, one = 1.0_dp
!
  integer               :: i, j
  logical               :: ok
  real(dp), allocatable :: eig(:), evec(:,:), evec_l(:,:)
  
!
! allocate space for eigenvalues and eigenvectors
!
  allocate (eig(n_max), evec(n, n_max))
  eig  = zero
  evec = zero
!
! test davidson:
!
  do i = 1, n_max
    evec(i,i) = one
  end do
  ok = .false.
!
  write(6,*) ' testing Davidson:'
  call davidson_driver(.true., n, n_targ, n_max, max_iter, tol, max_dav, shift, &
                       ax, dx, eig, evec, ok)
!
  if (ok) then
    write(6,*) ' Davidson converged.'
    write(6,*) ' Eigenvalues:'
    do i = 1, n_targ
      write(6,'(f16.8)') eig(i)
    end do
  else
    write(6,*) ' Davidson failed to converge.'
  end if
!
! test non-symmetric davidson:
!
  allocate (evec_l(n,n_max))
  evec_l = zero
  do i = 1, n_max
    evec_l(i,i) = one
  end do
  ok = .false.
!
  write(6,*) ' testing non-symmetric Davidson:'
  call nonsym_driver(.true., n, n_targ, n_max, max_iter, tol, max_dav, shift, &
                     arx, alx, dx, eig, evec, evec_l, 3, ok)
!
  if (ok) then
    write(6,*) ' non-symmetric Davidson converged.'
    write(6,*) ' Eigenvalues:'
    do i = 1, n_targ
      write(6,'(f16.8)') eig(i)
    end do
  else
    write(6,*) ' non-symmetric Davidson failed to converge.'
  end if
!
  deallocate (evec_l)
!
! test lobpcg:
!
  eig  = zero
  evec = zero
!
  do i = 1, n_max
    evec(i,i) = one
  end do
  ok = .false.
!
  write(6,*) ' testing LOBPCG:'
  call lobpcg_driver(.true., .true., n, n_targ, n_max, max_iter, tol, shift, &
                     ax, dx, sx, eig, evec, ok)
!
  if (ok) then
    write(6,*) ' LOBPCG converged.'
    write(6,*) ' Eigenvalues:'
    do i = 1, n_targ
      write(6,'(f16.8)') eig(i)
    end do
  else
    write(6,*) ' LOBPCG failed to converge.'
  end if
!
! test smogd:
!
  deallocate (evec)
  allocate (evec(2*n, n_max))
  eig  = zero
  evec = zero
!
  do i = 1, n_max
    evec(i,i)   = one
  end do
  ok = .false.
!
  write(6,*) ' testing SMOGD:'
  call smogd_driver(.true., n, 2*n, n_targ, n_max, max_iter, tol, max_dav, &
                    apbx, ambx, spdx, smdx, lrprc, eig, evec, ok)
  if (ok) then
    write(6,*) ' SMOGD converged.'
    write(6,*) ' Eigenvalues:'
    do i = 1, n_targ
      write(6,'(f16.8)') eig(i)
    end do
  else
    write(6,*) ' SMOGD failed to converge.'
  end if
!
! free the memory:
!
  deallocate (evec, eig)
!
  contains
!
  subroutine ax(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer :: i, j, k
!
    y = 0.0_dp
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + real(i+1,dp) * x(j,k)
          else
            y(i,k) = y(i,k) + x(j,k) / real(i+j,dp)
          end if
        end do
      end do
    end do
!
    return
  end subroutine ax
!
  subroutine dx(n,m,shift,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: shift
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, k
    real(dp) :: fac
!
    real(dp), parameter :: eps = 1.0e-5_dp
!
    do k = 1, m
      do i = 1, n
        fac = shift + real(i+1,dp)
        if (abs(fac).gt.eps) then
          y(i,k) = x(i,k) / (shift + real(i+1,dp))
        else
          y(i,k) = x(i,k)
        end if
      end do
    end do
  end subroutine dx
!
  subroutine arx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
    real(dp) :: fac
!
    y = 0.0_dp
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + real(i+1,dp) * x(j,k)
          else
            fac = real(i,dp) / real(j,dp)
            y(i,k) = y(i,k) + fac * x(j,k) / real(i+j,dp)
          end if
        end do
      end do
    end do
!
    return
  end subroutine arx
!
  subroutine alx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
    real(dp) :: fac
!
    y = 0.0_dp
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + real(i+1,dp) * x(j,k)
          else
            fac = real(j,dp) / real(i,dp)
            y(i,k) = y(i,k) + fac * x(j,k) / real(i+j,dp)
          end if
        end do
      end do
    end do
!
    return
  end subroutine alx
!
  subroutine sx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
!   just the identity matrix.
!
    y = x
!
    return
  end subroutine sx
!
  subroutine apbx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   (a + b)_ij = (5 + i) \delta_ij + (1 - \delta_ij) / (i+j)
!
    y = 0.0_dp
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + real(5+i,dp) * x(i,k)
          else
            y(i,k) = y(i,k) + x(j,k) / real(i+j)
          end if
        end do
      end do
    end do
    return        
  end subroutine apbx
!
  subroutine ambx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   (a + b)_ij = (2 + i) \delta_ij + (0.2 - \delta_ij) / (i+j)
!
    y = 0.0_dp
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + real(2+i,dp) * x(i,k)
          else
            y(i,k) = y(i,k) + 0.20_dp * x(j,k) / real(i+j)
          end if
        end do
      end do
    end do
    return        
  end subroutine ambx
!
  subroutine spdx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   \sigma = 1, \delta_ij = +- 0.05 
!
    y = 0.0_dp
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + x(i,k)
          else if (i.gt.j) then
            y(i,k) = y(i,k) + 0.05_dp * x(j,k) 
          else
            y(i,k) = y(i,k) - 0.05_dp * x(j,k) 
          end if
        end do
      end do
    end do
    return        
  end subroutine spdx
!
  subroutine smdx(n,m,x,y)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   (a + b)_ij = (2 + i) \delta_ij + (0.2 - \delta_ij) / (i+j)
!
    y = 0.0_dp
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + x(i,k)
          else if (i.gt.j) then
            y(i,k) = y(i,k) - 0.05_dp * x(j,k) 
          else
            y(i,k) = y(i,k) + 0.05_dp * x(j,k) 
          end if
        end do
      end do
    end do
    return        
  end subroutine smdx
!
  subroutine lrprc(n,m,fac,xp,xm,yp,ym)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer  :: i, k
    real(dp) :: val
!
!   yp = xp
!   ym = xm
    do k = 1, m
      do i = 1, n
        val = fac * fac * (real(i+7)**2 - 1.0_dp)
        val = 1.0_dp / val
        yp(i,k) = val * (fac * real(i+7) * xp(i,k) + xm(i,k))
        ym(i,k) = val * (fac * real(i+7) * xm(i,k) + xp(i,k))
      end do
    end do
!
    return
  end subroutine lrprc
end program test_fortran

