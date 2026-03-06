module direct_matvecs
use dgl_interface, only: dgl_real
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

end module direct_matvecs