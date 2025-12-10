program test_fortran
  use dgl_global_utils, only: dp, zero, one
  use dgl_drivers_interfaces
  implicit none
!
! tests the various functionalities of diaglib.
!
  integer,  parameter :: n = 100, n_targ = 5, n_max = 10, max_iter = 100, max_dav = 20
  integer,  parameter :: lutest=100
  real(dp), parameter :: tol = 1.0e-10_dp, shift = 0.0_dp
!
  integer               :: i, j, memory = 1.e7
  logical               :: ok
  real(dp), allocatable :: eig(:), evec(:,:), evec_l(:,:)
  procedure(), pointer :: mx_p => null()
!
! open a text file for the output, to be used to compare the results with a reference.
!
  open (file='output_fortran.txt',form='formatted',access='sequential', & 
        unit=lutest,status='unknown')
!
  1000 format(t3,a,1x,'results')
  1010 format(t3,'Eigenvalues:')
  1020 format(t3,i5,f14.6)
  1021 format(t3,i5,f12.4)
  1022 format(t3,i5,*(f14.6))
  1030 format(t3,'Eigenvector ',i3,':')
  1031 format(t3,'Eigenvector ',tl1,*(i8,6x))
!
! allocate space for eigenvalues and eigenvectors
!
  allocate (eig(n_max), evec(n, n_max))
!
! test davidson:
!
  eig  = zero
  evec = zero
  do i = 1, n_max
    evec(i,i) = one
  end do
  ok = .false.
!
  mx_p => mx
  write(6,*) ' testing Davidson:'
  call davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                        dgl_verbose = .true.)
!
  if (ok) then
    write(6,*) ' Davidson converged.'
    write(lutest,1000) 'Davidson'
    write(lutest,*)
    write(lutest,1010)
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do
    do i = 1, n_targ
!
!     fix the phase
!
      if (evec(1,i).lt.zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) ' Davidson failed to converge.'
  end if
!
! test generalized davidson:
!
  eig  = zero
  evec = zero
  do i = 1, n_max
    evec(i,i) = one
  end do
  ok = .false.
!
  write(6,*) ' testing Generalized Davidson:'
  call davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p)
!
  if (ok) then
    write(6,*) ' Generalized Davidson converged.'
    write(lutest,1000) 'Generalized Davidson'
    write(lutest,*)
    write(lutest,1010)
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do
    do i = 1, n_targ
!
!     fix the phase
!
      if (evec(1,i).lt.zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) 'Generalized Davidson failed to converge.'
  end if
!!
!! test non-symmetric davidson:
!!
!  allocate (evec_l(n,n_max))
!  evec_l = zero
!  do i = 1, n_max
!    evec_l(i,i) = one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing non-symmetric Davidson:'
!  call nonsym_driver(.false., n, n_targ, n_max, max_iter, tol, max_dav, shift, &
!                     arx, alx, dx, eig, evec, evec_l, 3, ok)
!!
!  if (ok) then
!    write(6,*) ' non-symmetric Davidson converged.'
!    write(lutest,1000) 'Non-Symmetric Davidson'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!    do i = 1, n_targ
!      write(lutest,1030) 'right', i
!!
!!     fix the phase
!!
!      if (evec(1,i).lt.zero) evec(:,i) = - evec(:,i)
!      do j = 1, n
!        write(lutest,1021) j, evec(j,i)
!      end do
!    end do
!    do i = 1, n_targ
!      write(lutest,1030) 'left', i
!!
!!     fix the phase
!!
!      if (evec_l(1,i).lt.zero) evec_l(:,i) = - evec_l(:,i)
!      do j = 1, n
!        write(lutest,1021) j, evec_l(j,i)
!      end do
!    end do
!    write(lutest,*)
!  else
!    write(6,*) ' non-symmetric Davidson failed to converge.'
!  end if
!!
!  deallocate (evec_l)
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

  call lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok)
!
  if (ok) then
    write(6,*) ' LOBPCG converged.'
    write(lutest,1000) 'LOBPCG'
    write(lutest,*)
    write(lutest,1010)
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do
    do i = 1, n_targ
!
!     fix the phase
!
      if (evec(1,i).lt.zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) ' LOBPCG failed to converge.'
  end if
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
  write(6,*) ' testing Generalized LOBPCG:'

  call lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p)
!
  if (ok) then
    write(6,*) ' Generalized LOBPCG converged.'
    write(lutest,1000) 'Generalized LOBPCG'
    write(lutest,*)
    write(lutest,1010)
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do
    do i = 1, n_targ
!
!     fix the phase
!
      if (evec(1,i).lt.zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) ' Generalized LOBPCG failed to converge.'
  end if
!!
!! test smogd:
!!
!  deallocate (evec)
!  allocate (evec(2*n, n_max))
!  eig  = zero
!  evec = zero
!!
!  do i = 1, n_max
!    evec(i,i)   = one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing SMOGD:'
!  call smogd_driver(.false., n, 2*n, n_targ, n_max, max_iter, tol, max_dav, &
!                    apbx, ambx, spdx, smdx, lrprc, eig, evec, ok)
!  if (ok) then
!    write(6,*) ' SMOGD converged.'
!    write(lutest,1000) 'SMOGD'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!    write(lutest,*)
!  else
!    write(6,*) ' SMOGD failed to converge.'
!  end if
!
! close the output file:
!
  close (lutest)
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
  subroutine mx(n,m,x,y)
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
            y(i,k) = y(i,k) + x(j,k)
          else
            y(i,k) = y(i,k) + x(j,k) / real(i+j,dp)
          end if
        end do
      end do
    end do
!
    return
  end subroutine mx
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

  subroutine test(a)
    real :: a
    print *, "is called, a=", a
  end subroutine
end program test_fortran

