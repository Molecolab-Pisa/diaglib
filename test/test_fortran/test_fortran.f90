program test_fortran
  use dgl_interface
  implicit none
!
! tests the various functionalities of diaglib.
!
  integer,  parameter :: n = 500, n_targ = 5, n_max = 10, max_iter = 100, max_dav = 20
  integer,  parameter :: lutest=100
  real(dgl_real), parameter :: tol = 1.0e-10_dgl_real, shift = 0.0_dgl_real
  integer,  parameter :: memory = 100
  character(len=2),parameter :: memory_unit = "MB"
!
  integer               :: i, j
  logical               :: ok
  real(dgl_real), allocatable :: eig(:), evec(:,:), evec_l(:,:)
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
  1032 format(t3,a6,'Eigenvector ',tl7,*(i8,6x))
!
! allocate space for eigenvalues and eigenvectors
!
  allocate (eig(n_max), evec(n, n_max))
!
! test davidson:
!
  eig  = dgl_zero
  evec = dgl_zero
  do i = 1, n_max
    evec(i,i) = dgl_one
  end do
  ok = .false.
!
  mx_p => mx
  write(6,*) ' testing Davidson:'
  call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                        dgl_verbose = .false.,&
                        dgl_memory = memory,&
                        dgl_memory_unit = memory_unit)
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
  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
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
  eig  = dgl_zero
  evec = dgl_zero
  do i = 1, n_max
    evec(i,i) = dgl_one
  end do
  ok = .false.
!
  write(6,*) ' testing Generalized Davidson:'
  call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
                        dgl_verbose = .false.,&
                        dgl_memory = memory,&
                        dgl_memory_unit = memory_unit)
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
! fix the phase
!
  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) 'Generalized Davidson failed to converge.'
  end if
!
! test non-symmetric davidson:
!
  allocate (evec_l(n,n_max))
  evec = dgl_zero
  do i = 1, n_max
    evec(i,i) = dgl_one
  end do
  ok = .false.
!
  write(6,*) ' testing non-symmetric Davidson:'
  call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "LR", eig, evec, ok, evec_2 = evec_l, &
                              dgl_verbose = .false.,&
                              dgl_memory = memory, &
                              dgl_memory_unit = memory_unit)
!
  if (ok) then
    write(6,*) ' non-symmetric Davidson converged.'
    write(lutest,1000) 'Non-Symmetric Davidson'
    write(lutest,*)
    write(lutest,1010)
    
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do

    write(lutest,1032) 'Right ', (i, i = 1, n_targ)
!
!   fix the phase
!
    if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)

    write(lutest,1032) 'Left ', (i, i = 1, n_targ)
!
!     fix the phase
!
    if (evec_l(1,i).lt.dgl_zero) evec_l(:,i) = - evec_l(:,i)
    do j = 1, n
      write(lutest,1022) j, (evec_l(j,i), i = 1, n_targ)
    end do
    write(lutest,*)

  else
    write(6,*) ' non-symmetric Davidson failed to converge.'
  end if
!
  deallocate (evec_l)
!
! test lobpcg:
!
  eig  = dgl_zero
  evec = dgl_zero
!
  do i = 1, n_max
    evec(i,i) = dgl_one
  end do
  ok = .false.
!
  write(6,*) ' testing LOBPCG:'

  call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                        dgl_verbose = .false.,&
                        dgl_memory = memory,&
                        dgl_memory_unit = memory_unit)
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
  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
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
  eig  = dgl_zero
  evec = dgl_zero
!
  do i = 1, n_max
    evec(i,i) = dgl_one
  end do
  ok = .false.
!
  write(6,*) ' testing Generalized LOBPCG:'

  call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
                        dgl_verbose = .false.,&
                        dgl_memory = memory,&
                        dgl_memory_unit = memory_unit)
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
  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) ' Generalized LOBPCG failed to converge.'
  end if
!
! test smogd:
!
  deallocate (evec)
  allocate (evec(2*n, n_max))
  eig  = dgl_zero
  evec = dgl_zero
!
  do i = 1, n_max
    evec(i,i)   = dgl_one
  end do
  ok = .false.
!
  write(6,*) ' testing SMOGD:'
  call dgl_smogd_driver(2*n, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, &
                        eig, evec, ok, &
                        dgl_dav_iter = 5, &
                        dgl_tol = 1.e-12_dgl_real, &
                        dgl_max_iter = 20, &
                        dgl_verbose = .true., &
                        dgl_memory = memory, &
                        dgl_memory_unit = memory_unit)
  if (ok) then
    write(6,*) ' SMOGD converged.'
    write(lutest,1000) 'SMOGD'
    write(lutest,*)
    write(lutest,1010)
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do
    write(lutest,*)
  else
    write(6,*) ' SMOGD failed to converge.'
  end if
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer :: i, j, k
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + real(i+1,dgl_real) * x(j,k)
          else
            y(i,k) = y(i,k) + x(j,k) / real(i+j,dgl_real)
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer :: i, j, k
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + x(j,k)
          else
            y(i,k) = y(i,k) + x(j,k) / real(i+j,dgl_real)
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
    real(dgl_real),                 intent(in)    :: shift
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, k
    real(dgl_real) :: fac
!
    real(dgl_real), parameter :: eps = 1.0e-5_dgl_real
!
    do k = 1, m
      do i = 1, n
        fac = shift + real(i+1,dgl_real)
        if (abs(fac).gt.eps) then
          y(i,k) = x(i,k) / (shift + real(i+1,dgl_real))
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
    real(dgl_real) :: fac
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + real(i+1,dgl_real) * x(j,k)
          else
            fac = real(i,dgl_real) / real(j,dgl_real)
            y(i,k) = y(i,k) + fac * x(j,k) / real(i+j,dgl_real)
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
    real(dgl_real) :: fac
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do i = 1, n
        do j = 1, n
          if (j.eq.i) then
            y(i,k) = y(i,k) + real(i+1,dgl_real) * x(j,k)
          else
            fac = real(j,dgl_real) / real(i,dgl_real)
            y(i,k) = y(i,k) + fac * x(j,k) / real(i+j,dgl_real)
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   (a + b)_ij = (5 + i) \delta_ij + (1 - \delta_ij) / (i+j)
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + real(5+i,dgl_real) * x(i,k)
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   (a + b)_ij = (2 + i) \delta_ij + (0.2 - \delta_ij) / (i+j)
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + real(2+i,dgl_real) * x(i,k)
          else
            y(i,k) = y(i,k) + 0.20_dgl_real * x(j,k) / real(i+j)
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   \sigma = 1, \delta_ij = +- 0.05 
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + x(i,k)
          else if (i.gt.j) then
            y(i,k) = y(i,k) + 0.05_dgl_real * x(j,k) 
          else
            y(i,k) = y(i,k) - 0.05_dgl_real * x(j,k) 
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
    real(dgl_real), dimension(n,m), intent(in)    :: x
    real(dgl_real), dimension(n,m), intent(inout) :: y
!
    integer  :: i, j, k
!
!   (a + b)_ij = (2 + i) \delta_ij + (0.2 - \delta_ij) / (i+j)
!
    y = 0.0_dgl_real
!
    do k = 1, m
      do j = 1, n
        do i = 1, n
          if (i.eq.j) then
            y(i,k) = y(i,k) + x(i,k)
          else if (i.gt.j) then
            y(i,k) = y(i,k) - 0.05_dgl_real * x(j,k) 
          else
            y(i,k) = y(i,k) + 0.05_dgl_real * x(j,k) 
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
    real(dgl_real),                 intent(in)    :: fac
    real(dgl_real), dimension(n,m), intent(in)    :: xp, xm
    real(dgl_real), dimension(n,m), intent(inout) :: yp, ym
!
    integer  :: i, k
    real(dgl_real) :: val
!
!   yp = xp
!   ym = xm
    do k = 1, m
      do i = 1, n
        val = fac * fac * (real(i+7)**2 - 1.0_dgl_real)
        val = 1.0_dgl_real / val
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

