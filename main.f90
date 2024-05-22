program main
  use real_precision
  use utils
  use diaglib
!
! simple minded test program to call the lobpcg routine
!
  implicit none
  integer  :: n, itmax, m_max, n_want, iwhat, cwhat, pwhat
  real(dp) :: tol
!
! initialize:
!
  n      = 200
  n_want = 10
  tol    = 1.0e-8_dp
  itmax  = 1000
  m_max  = 10
  nmult  = 0
  tdscf  = .false.
  i_alg  = 0
!
  if ((m_max*n_want) .gt. n) then 
     write(6,*) 
     write(6,*) 'WARNING' 
     write(6,*) 'Rouché-Capelli is encountered! Please, change your parameters.' 
     write(6,*) 
     stop
   end if
!
  call print_header
! 
  write(6,*)
  1000 format(t3,' simple-minded test driver. select ',/, &
              t3,'   1 for symmetric eigenvalue problems ',/, &
              t3,'   2 for symmetric generalized eigenvalue problems ',/, &
              t3,'   3 for linear-response equations (SCF-like)',/, &
              t3,'   4 for linear-response equations (CASSCF-like)',/, &
              t3,'   5 for linear-response equations (complex-CASSCF-like)',/, &
              t3,'   6 for symmetric eigenvalue problems (complex-Davidson).')
  write(6,1000)
  read(5,*) iwhat
  !iwhat = 5
  write(6,*)
!
  if (iwhat.eq.1) then 
    call test_symm(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.2) then 
    call test_geneig(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.3) then 
    call test_scflr(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.4) then 
    call test_caslr(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.5) then 
    call test_caslr_newcomplex(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.6) then 
    call test_symm_complex(.true.,n,n_want,tol,itmax,m_max)
  else 
    write(6,*) ' invalid selectrion. aborting ...'
  end if
!
end program main
!
  subroutine print_header
    implicit none
!
    1000 format( &
  ' diaglib - a fortran library of matrix-free iterative algorithms to ',/, &
  ' compute a few eigenvalues and eigenvectors of large matrices.',/, &
  ' ==================================================================',/, &
  /, &
  ' Implementation by',/, &
  /, &
  '   Ivan Giannì, Tommaso Nottoli, Federica Pes, Antoine Levitt ',/, &
  '   and Filippo Lipparini',/, &
  '   MoLECoLab Pisa',/, &
  '   Department of Chemistry and Industrial Chemistry',/, &
  '   University of Pisa',/, &
  '   Via G. Moruzzi 13, I-56124, Pisa, Italy')
    write(6,1000)
    return
  end subroutine print_header
!
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
  subroutine mmult_newcomplex(n,m,x,ax)
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
    integer :: icol, reim
!
    nmult = nmult + m
    do icol = 1, m
      ax(:,icol) = matmul(fourmat,x(:,icol))
    end do
    return
  end subroutine mmult_newcomplex
!
  subroutine mmult_complex(n,m,x,ax)
    use utils
    implicit none
!
!   simple-minded matrix-vector multiplication routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                      intent(in)    :: n, m
    complex(dp),  dimension(n,m), intent(in)    :: x
    complex(dp),  dimension(n,m), intent(inout) :: ax
!
    integer :: icol
!
    nmult = nmult + m
    do icol = 1, m
      ax(:,icol) = matmul(ca,x(:,icol))
    end do
    return
  end subroutine mmult_complex
!
 subroutine smult(n,m,x,sx)
    use utils
    implicit none
!
!   simple-minded matrix-vector multiplication routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,m), intent(inout) :: sx
!
    integer :: icol
!
!
!
    if (tdscf) then
!
!     nothing to do, as s is the identity. not the most elegant way of doing this,
!     but it works.
!
      sx = x
      return
    end if
!
    nmult = nmult + m
    do icol = 1, m
      sx(:,icol) = matmul(s,x(:,icol))
    end do
    return
  end subroutine smult
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
  subroutine mprec_newcomplex(n,m,fac,x,px)
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
        if (abs(fourmat(i,i)+fac).gt.tol) px(i,icol) = x(i,icol)/(fourmat(i,i) + fac)
      end do
    end do
    return
  end subroutine mprec_newcomplex
!

  subroutine mprec_complex(n,m,fac,x,px)
    use utils
    implicit none
!
!   simple-minded shift-and-invert preconditioner routine, that needs to be passed as 
!   an argument to lobpcg
!
    integer,                      intent(in)    :: n, m
    complex(dp),                  intent(in)    :: fac
    complex(dp),  dimension(n,m), intent(in)    :: x
    complex(dp),  dimension(n,m), intent(inout) :: px
!
    integer             :: i, icol
    real(dp), parameter :: tol = 1.0d-5
!
    do icol = 1, m
      do i = 1, n
        if (abs(ca(i,i)+fac).gt.tol) px(i,icol) = x(i,icol)/(ca(i,i) + fac)
      end do
    end do
    return
  end subroutine mprec_complex
!
  subroutine lambdavec(n,m,x,y)
    use utils
    implicit none 
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(lambdafull,x(:,icol))
    end do
!
    return
  end subroutine lambdavec
!
  subroutine omegavec(n,m,x,y)
    use utils
    implicit none 
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(omegafull,x(:,icol))
    end do
!
    return
  end subroutine omegavec
!
  subroutine apbvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(apb,x(:,icol))
    end do
!
    return
  end subroutine apbvec
!
  subroutine ambvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(amb,x(:,icol))
    end do
!
    return
  end subroutine ambvec
!
  subroutine spdvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(spd,x(:,icol))
    end do
!
    return
  end subroutine spdvec
!
  subroutine smdvec(n,m,x,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x
    real(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
    do icol = 1, m
      y(:,icol) = matmul(smd,x(:,icol))
    end do
!
    return
  end subroutine smdvec
!
  subroutine apbvec_newcomplex(reim,n,m,x1,x2,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x1, x2
    real(dp), dimension(n,m), intent(inout) :: y
!
!   if reim = 0, use apbre; if reim = 1 use apbim
!
    integer                                 :: reim
!
    integer icol
!
    if (reim.eq.0) then
!
!     here we want to compute sigma_re+ = (A+B)re*Vre+ - (A-B)im*Vim+
!
!      write(6,*) 'hereeee apbre'
      do icol = 1, m
        y(:,icol) = matmul(apbre,x1(:,icol)) - matmul(ambim,x2(:,icol))
      end do
    else if (reim.eq.1) then
!
!     here we want to compute sigma_im+ = (A+B)im*Vre+ - (A-B)re*Vim+
!
!      write(6,*) 'hereeee apbim'
      do icol = 1, m
        y(:,icol) = matmul(apbim,x1(:,icol)) + matmul(ambre,x2(:,icol))
      end do
    else
      write(6,*) 'invalid entry: first position of apbvec_newcomplex subroutine must be either 0 or 1'
    end if
!
    return
  end subroutine apbvec_newcomplex
!
  subroutine ambvec_newcomplex(reim,n,m,x1,x2,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x1,x2
    real(dp), dimension(n,m), intent(inout) :: y
!
!   if reim = 0, use ambre; if reim = 1 use ambim
!
    integer                                 :: reim
!
    integer icol
!
    if (reim.eq.0) then
!
!     here we want to compute sigma_re- = (A-B)re*Vre- - (A+B)im*Vim-
!
!      write(6,*) 'hereeee ambre'
      do icol = 1, m
        y(:,icol) = matmul(ambre,x1(:,icol)) - matmul(apbim,x2(:,icol))
      end do
    else if (reim.eq.1) then
!
!     here we want to compute sigma_im- = (A-B)im*Vre- + (A+B)re*Vim-
!
!    write(6,*) 'hereeee ambim'
      do icol = 1, m
        y(:,icol) = matmul(ambim,x1(:,icol)) + matmul(apbre,x2(:,icol))
      end do
    else
      write(6,*) 'invalid entry: first position of ambvec_newcomplex subroutine must be either 0 or 1'
    end if
!
    return
  end subroutine ambvec_newcomplex
!
  subroutine spdvec_newcomplex(reim,n,m,x1,x2,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x1, x2
    real(dp), dimension(n,m), intent(inout) :: y
!
!   if reim = 0, use spdre; if reim = 1 use spdim
!
    integer                                 :: reim
!
    integer icol
!
    if (reim.eq.0) then
!
!     here we want to compute tau_re- = (SIG + DEL)re*Vre+ - (SIG - DEL)im*Vim+
!
!      write(6,*) 'hereeee spdre'
      do icol = 1, m
        y(:,icol) = matmul(spdre,x1(:,icol)) - matmul(smdim,x2(:,icol))
      end do
      !write(6,*) 'bvmre in sub', y
    else if (reim.eq.1) then
!
!     here we want to compute tau_im- = (SIG + DEL)im*Vre+ + (SIG - DEL)re*Vim+
!
!      write(6,*) 'hereeee spdim'
      do icol = 1, m
        y(:,icol) = matmul(spdim,x1(:,icol)) + matmul(smdre,x2(:,icol))
      end do
    else
      write(6,*) 'invalid entry: first position of spdvec_newcomplex subroutine must be either 0 or 1'
    end if
!
    return
  end subroutine spdvec_newcomplex
!
  subroutine smdvec_newcomplex(reim,n,m,x1,x2,y)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(in)    :: x1, x2
    real(dp), dimension(n,m), intent(inout) :: y
!
!   if reim = 0, use smdre; if reim = 1 use smdim
!
    integer                                 :: reim
!
    integer icol
!
    if (reim.eq.0) then
!
!     here we want to compute tau_re+ = (SIG - DEL)re*Vre- - (SIG + DEL)im*Vim-
!
!      write(6,*) 'hereeee smdre'
      do icol = 1, m
        y(:,icol) = matmul(smdre,x1(:,icol)) - matmul(spdim,x2(:,icol))
      end do
    else if (reim.eq.1) then
!
!     here we want to compute tau_im+ = (SIG - DEL)im*Vre- + (SIG + DEL)re*Vim-
!
!      write(6,*) 'hereeee smdim'
      do icol = 1, m
        y(:,icol) = matmul(smdim,x1(:,icol)) + matmul(spdre,x2(:,icol))
      end do
    else
      write(6,*) 'invalid entry: first position of smdvec_newcomplex subroutine must be either 0 or 1'
    end if
!
    return
  end subroutine smdvec_newcomplex
!
  subroutine apbvec_complex(n,m,x,y)
    use utils
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp), dimension(n,m), intent(in)    :: x
    complex(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
!   we want to compute sigma_i^+ =Ab^+_i + B(b^+_i)* 
!
    do icol = 1, m
      y(:,icol) = matmul(caa,x(:,icol)) + matmul(cbb,conjg(x(:,icol)))
    end do
!
    return
  end subroutine apbvec_complex
!
  subroutine ambvec_complex(n,m,x,y)
    use utils
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp), dimension(n,m), intent(in)    :: x
    complex(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
!   we want to compute sigma_i^- =Ab^-_i - B(b^-_i)* 
!
!    write(6,*) 'x_amb', x
    do icol = 1, m
      y(:,icol) = matmul(caa,x(:,icol)) - matmul(cbb,conjg(x(:,icol)))
    end do
!
    return
  end subroutine ambvec_complex
!
  subroutine spdvec_complex(n,m,x,y)
    use utils
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp), dimension(n,m), intent(in)    :: x
    complex(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
!   we want to compute tau_i^- = Sigmab^+_i + Delta(b^+_i)* 
!
!    write(6,*) 'csigma', csigma 
!    write(6,*) 'cdelta', cdelta 
    do icol = 1, m
      y(:,icol) = matmul(csigma,x(:,icol)) + matmul(cdelta,conjg(x(:,icol)))
    end do
    !call zgemm('n','n',n,m,n,(1.0_dp,0.0_dp),csigma,n,x,n,(0.0_dp,0.0_dp),y,n)
    !call zgemm('n','n',n,m,n,(1.0_dp,0.0_dp),cdelta,n,conjg(x),n,(1.0_dp,0.0_dp),y,n)
!
    return
  end subroutine spdvec_complex
!
  subroutine smdvec_complex(n,m,x,y)
    use utils
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp), dimension(n,m), intent(in)    :: x
    complex(dp), dimension(n,m), intent(inout) :: y
!
    integer icol
!
!   we want to compute tau_i^+ = Sigmab^-_i - Delta(b^-_i)* 
!
!    write(6,*) 'caa', caa
!    write(6,*) 'cdelta', cdelta
!    write(6,*) 'csigma', csigma
!    write(6,*) 'x_vm', x
    do icol = 1, m
      y(:,icol) = matmul(csigma,x(:,icol)) - matmul(cdelta,conjg(x(:,icol)))
    end do
    !call zgemm('n','n',n,m,n,(1.0_dp,0.0_dp),csigma,n,x,n,(0.0_dp,0.0_dp),y,n)
    !call zgemm('n','n',n,m,n,(1.0_dp,0.0_dp),cdelta,n,conjg(x),n,(-1.0_dp,0.0_dp),y,n)
!    write(6,*) 'ysmd', y
!
    return
  end subroutine smdvec_complex
!
  subroutine lrprec_2f(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer :: i, icol
!
    yp = 0.0_dp
    ym = 0.0_dp
!
    do icol = 1, m
      do i = 1, n
        yp(i,icol) = - 1.0_dp/(lambdafull(i,i)**2 - fac**2 * omegafull(i,i)**2) * & 
                              (lambdafull(i,i) * xp(i,icol) + fac * omegafull(i,i) * xm(i,icol))
        ym(i,icol) = - 1.0_dp/(lambdafull(i,i)**2 - fac**2 * omegafull(i,i)**2) * & 
                              (lambdafull(i,i) * xm(i,icol) + fac * omegafull(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_2f
!
  subroutine lrprec_1(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer :: i, icol
!
    yp = 0.0_dp
    ym = 0.0_dp
!
    do icol = 1, m
      do i = 1, n
        yp(i,icol) = - 1.0_dp/(aa(i,i)**2 - fac**2 * sigma(i,i)**2) * (aa(i,i) * xp(i,icol) + fac * sigma(i,i) * xm(i,icol))
        ym(i,icol) = - 1.0_dp/(aa(i,i)**2 - fac**2 * sigma(i,i)**2) * (aa(i,i) * xm(i,icol) + fac * sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_1
!
  subroutine lrprec_2(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer  :: i, icol
    real(dp) :: denom
!
    yp = 0.0_dp
    ym = 0.0_dp
!
    do icol = 1, m
      do i = 1, n
        denom = fac * fac * aa(i,i) * aa(i,i) - sigma(i,i) * sigma(i,i)
        denom = 1.0_dp/denom
        yp(i,icol) = denom * (fac * aa(i,i) * xp(i,icol) + sigma(i,i) * xm(i,icol))
        ym(i,icol) = denom * (fac * aa(i,i) * xm(i,icol) + sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_2
!
  subroutine lrprec_2_newcomplex(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer  :: i, icol
    real(dp) :: denom
!
    yp = 0.0_dp
    ym = 0.0_dp
!
    do icol = 1, m
      do i = 1, n
        denom = fac * fac * aare(i,i) * aare(i,i) - sigmare(i,i) * sigmare(i,i)
        denom = 1.0_dp/denom
        yp(i,icol) = denom * (fac * aare(i,i) * xp(i,icol) + sigmare(i,i) * xm(i,icol))
        ym(i,icol) = denom * (fac * aare(i,i) * xm(i,icol) + sigmare(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_2_newcomplex
!
  subroutine lrprec_1_newcomplex(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp),                 intent(in)    :: fac
    real(dp), dimension(n,m), intent(in)    :: xp, xm
    real(dp), dimension(n,m), intent(inout) :: yp, ym
    real(dp)                                :: denom, denom1
!
    integer :: i, icol
!
    yp = 0.0_dp
    ym = 0.0_dp
!
    do icol = 1, m
      do i = 1, n
        denom1     = (aare(i,i)**2 - fac**2 * sigmare(i,i)**2)
        denom      = 1.0_dp/denom1
        yp(i,icol) = - denom * (aare(i,i) * xp(i,icol) + fac * sigmare(i,i) * xm(i,icol))
        ym(i,icol) = - denom * (aare(i,i) * xm(i,icol) + fac * sigmare(i,i) * xp(i,icol))
      end do
    end do
    return
  end subroutine lrprec_1_newcomplex
!
  subroutine lrprec_1_complex(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp),                 intent(in)    :: fac
    complex(dp), dimension(n,m), intent(in)    :: xp, xm
    complex(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer     :: i, icol
    complex(dp) :: cdenom
!
    yp = (0.0_dp,0.0_dp)
    ym = (0.0_dp,0.0_dp)
!
    do icol = 1, m
      !print '(4a20)', 'Re(yp)', 'Im(yp)', 'Re(ym)', 'Im(ym)'
      do i = 1, n
        !yp(i,icol) = - (1.0_dp,0.0_dp)/(caa(i,i)**2 - fac*conjg(fac) * csigma(i,i)**2 * & 
        !               (caa(i,i) * xp(i,icol) + fac * csigma(i,i) * xm(i,icol)))
        !ym(i,icol) = - (1.0_dp,0.0_dp)/(caa(i,i)**2 - fac*conjg(fac) * csigma(i,i)**2 * &
        !               (caa(i,i) * xm(i,icol) + fac * csigma(i,i) * xp(i,icol)))
        cdenom     = caa(i,i)**2 - fac**2 * csigma(i,i)**2
        cdenom     = -(1.0_dp,0.0_dp)/cdenom
        yp(i,icol) = cdenom * (caa(i,i) * xp(i,icol) + fac * csigma(i,i) * xm(i,icol))
        ym(i,icol) = cdenom * (caa(i,i) * xm(i,icol) + fac * csigma(i,i) * xp(i,icol))
 !      ! print '(4f20.10)', real(yp(i,icol)), aimag(yp(i,icol)), real(ym(i,icol)), aimag(ym(i,icol))
!
! activate to test without precond  
!
!        yp(i,icol) = - xp(i,icol) 
!        ym(i,icol) = - xm(i,icol) 
      end do
      !stop
    end do

!
    return
  end subroutine lrprec_1_complex
!
  subroutine lrprec_2_complex(n,m,fac,xp,xm,yp,ym)
    use utils
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp),                 intent(in)    :: fac
    complex(dp), dimension(n,m), intent(in)    :: xp, xm
    complex(dp), dimension(n,m), intent(inout) :: yp, ym
!
    integer     :: i, icol
    real(dp)    :: denom
    complex(dp) :: cdenom
!
    yp = (0.0_dp,0.0_dp)
    ym = (0.0_dp,0.0_dp)
! 
    do icol = 1, m
      do i = 1, n
        cdenom = fac * fac * caa(i,i) * caa(i,i) - csigma(i,i) * csigma(i,i)
        cdenom = (1.0_dp,0.0_dp)/cdenom
        yp(i,icol) = cdenom * (fac * caa(i,i) * xp(i,icol) + csigma(i,i) * xm(i,icol))
        ym(i,icol) = cdenom * (fac * caa(i,i) * xm(i,icol) + csigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_2_complex
!
  subroutine test_symm(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : lobpcg_driver, davidson_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   test iterative solvers for a standard symmetric eigenvalue problem.
!
    logical  :: ok
    integer  :: i, j, n_eig, info, lwork
    real(dp) :: lw(1)
    real(dp), allocatable :: diagonal(:), eig(:), evec(:,:)
    real(dp), allocatable :: work(:), w(:), a_copy(:,:)
!
    external :: mmult, mprec, smult
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!
!   allocate memory for the a matrix 
!
    allocate (a(n,n))
!
!   build a symmetric matrix:
!
    do i = 1, n
      a(i,i) = real(i,kind=dp) + 1.0_dp
      do j = 1, i - 1
        a(j,i) = 1.0_dp/real(i+j,kind=dp)
        a(i,j) = a(j,i)
      end do
    end do
!
!   if required, solve the problem with a dense lapack routine:
!
    if (check_lapack) then
      allocate (a_copy(n,n), w(n))
      a_copy = a
      call dsyev('v','l',n,a_copy,n,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsyev('v','l',n,a_copy,n,w,work,lwork,info)
!
!     write the results on file for comparison:
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, w(i)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (a_copy(1,i) .lt. 0.0_dp) a_copy(:,i) = - a_copy(:,i)
        write(10,'(10f12.6)') a_copy(:,i)
        write(10,*)
      end do
      close (10)
    end if
!
!   allocate and gather the diagonal:
!
    allocate (diagonal(n))
    do i = 1, n
      diagonal(i) = a(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvalues and eigenvectors:
!
    allocate (eig(n), evec(n,n_eig))
!
!   compute a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call lobpcg:
!
    call lobpcg_driver(.true.,.false.,n,n_want,n_eig,itmax,tol,0.0_dp,mmult, &
                       mprec,smult,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'lobcpg.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call the davidson driver:
!
    call davidson_driver(.true.,n,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult, &
                         mprec,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'davidson.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_symm
!
  subroutine test_symm_complex(check_lapack,n,n_want,tol,itmax,m_max)
!          
!   complex davidson and new-complex davidson 
!
    use real_precision
    use utils
    use diaglib!, only : lobpcg_driver, davidson_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   test iterative solvers for a standard symmetric eigenvalue problem.
!
    logical                   :: ok
    integer                   :: i, j, n_eig, info, lwork, lrwork, liwork
    integer,     allocatable  :: iwork(:)
    real(dp)                  :: lw(1)
    real(dp),    allocatable  :: diagonal(:), eig(:), w(:), rwork(:) 
    complex(dp), allocatable  :: cwork(:), a_copy(:,:), ceig(:)
    complex(dp), allocatable  :: evec(:,:), cprod(:,:)
!
!   new complex scheme
!
    integer                   :: n2
    real(dp),    allocatable  :: nc_copy(:,:)
    real(dp),    allocatable  :: work(:), ncw(:)
    real(dp),    allocatable  :: ncdiag(:), nceig(:), evecre(:,:), evecim(:,:), ncevec(:,:)
!
    external :: mmult_complex, mprec_complex, smult, mmult_newcomplex, mprec_newcomplex
!
    write(6,*)'test davidson_complex'
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
    n2 = 2*n
!
!   allocate memory for the a matrix (ca is complex, are/im are for new-complex
!
    allocate (a(n,n), ca(n,n), are(n,n), aim(n,n), fourmat(n2,n2))
!
!   build a hermitian matrix:
!
    do i = 1, n
      ca(i,i)  = dcmplx(i+1.0_dp,0.0_dp) 
      are(i,i) = i + 1.0_dp
      do j = 1, i - 1
        ca(j,i)  = dcmplx(1.0_dp/(i+j),1.0_dp/(i+j))
        ca(i,j)  = conjg(ca(j,i))
        are(j,i) = (1.0_dp)/(i+j) 
        are(i,j) = are(j,i)
        aim(j,i) = (1.0_dp)/(i+j)
        aim(i,j) = -aim(j,i)
      end do
    end do
!
!   build the newcomplex matrix A=(Are -Aim, Aim Are)
!
    fourmat = 0.0_dp
    fourmat(1:n,1:n)       = are
    fourmat(n+1:n2,1:n)    = aim
    fourmat(1:n,n+1:n2)    = -aim
    fourmat(n+1:n2,n+1:n2) = are
!
!   if required, solve the problem with a dense lapack routine:
!
    if (check_lapack) then
      allocate (a_copy(n,n), w(n), cwork(1), rwork(1), iwork(1))
      a_copy = ca
      call zheevd('v','l',n,a_copy,n,w,cwork,-1,rwork,-1,iwork,-1,info)
      write(6,*) 'info1 cmplx', info
      lwork = int(cwork(1))
      lrwork = int(rwork(1))
      liwork = int(iwork(1))
      deallocate(cwork, rwork, iwork)
      allocate (cwork(lwork), rwork(lrwork), iwork(liwork))
      call zheevd('v','l',n,a_copy,n,w,cwork,lwork,rwork,lrwork,iwork,liwork,info)
      write(6,*) 'info2 cmplx', info
!
!     solve the newcomplex problem with a dense lapack routine:
!
      allocate(nc_copy(n2,n2), ncw(n2))
      nc_copy = fourmat
      call dsyev('v','l',n2,nc_copy,n2,ncw,lw,-1,info)
      write(6,*) 'info1 newcomplex', info
      lwork = int(lw(1))
      allocate(work(lwork))
      call dsyev('v','l',n2,nc_copy,n2,ncw,work,lwork,info)
      write(6,*) 'info2 newcomplex', info
!
!     write the results on file for comparison:
!
      open (unit = 10, file = 'lapack_dav_complex.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, w(i)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (abs(a_copy(1,i)) .lt. 0.0_dp) a_copy(:,i) = - a_copy(:,i)
        write(10,*) a_copy(:,i)
        write(10,*)
      end do
      close (10)
!
      open (unit = 10, file = 'lapack_dav_newcomplex.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, ncw(i)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (abs(nc_copy(1,i)) .lt. 0.0_dp) nc_copy(:,i) = - nc_copy(:,i)
        write(10,*) nc_copy(:,i)
        write(10,*)
      end do
      close (10)
    end if
!
!   allocate and gather the diagonal:
!
    allocate (diagonal(n), ncdiag(n2))
    do i = 1, n
      diagonal(i) = real(ca(i,i))
      ncdiag(i)   = fourmat(i,i)
      ncdiag(i+n) = fourmat(i+n,i+n)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
    n_eig = n_want
!
!   allocate memory for the eigenvalues and eigenvectors:
!
    allocate (eig(n), ceig(n), evec(n,n_eig), cprod(n,n_eig))
    allocate (nceig(n2), ncevec(n2,n_eig), evecre(n,n_eig), evecim(n,n_eig))
!
!   complex davidson
!
    write(6,*) 
    write(6,*) '----------------------' 
    write(6,*) '   COMPLEX DAVIDSON'
    write(6,*) '----------------------' 
    write(6,*) 
!
!   make a guess for the eigenvector
!
    call guess_evec_complex(6,n,n_eig,diagonal,evec)
!
!   call the davidson driver:
!
    call davidson_complex_driver(.true.,n,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult_complex, &
                                 mprec_complex,eig,ceig,evec,ok)
!
!   test A x=lambda x
!
    !write(6,*) 'TEST AX=LAMBDAX'
    !call zgemv('n',n,n,(1.0_dp,0.0_dp),ca,n,evec(:,1),1,(0.0_dp,0.0_dp),cprod,1)
    !write(6,*) 'cprod', cprod
    !write(6,*) 'lambdax', eig(1)*evec(:,1)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'davidson_complex.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (abs(evec(1,i)) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,*) evec(:,i)
      write(10,*)
    end do
    close (10)
!
!   new complex davidson
!
    write(6,*) 
    write(6,*) '--------------------------' 
    write(6,*) '   NEW COMPLEX DAVIDSON'
    write(6,*) '--------------------------' 
    write(6,*) 
!
!   make a guess for the eigenvector
!
    call guess_evec(4,n2,n_eig,ncdiag,ncevec)
!
!   call the davidson driver:
!
    !call davidson_newcomplex_driver(.true.,n2,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult_newcomplex, &
    !                             mprec_newcomplex,nceig,ncevec,ok)
    call davidson_driver(.true.,n2,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult_newcomplex, &
                         mprec_newcomplex,nceig,ncevec,ok)
!
!   test A x=lambda x
!
    !write(6,*) 'TEST AX=LAMBDAX'
    !call zgemv('n',n,n,(1.0_dp,0.0_dp),ca,n,evec(:,1),1,(0.0_dp,0.0_dp),cprod,1)
    !write(6,*) 'cprod', cprod
    !write(6,*) 'lambdax', eig(1)*evec(:,1)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'davidson_newcomplex.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, nceig(i)
      if (abs(ncevec(1,i)) .lt. 0.0_dp) ncevec(:,i) = - ncevec(:,i)
      write(10,*) ncevec(:,i)
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_symm_complex
!
  subroutine test_geneig(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : lobpcg_driver, gen_david_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
    logical  :: ok
    integer  :: i, j, n_eig, info, lwork
    real(dp) :: lw(1)
    real(dp), allocatable :: diagonal(:), eig(:), evec(:,:)
    real(dp), allocatable :: work(:), w(:), a_copy(:,:), s_copy(:,:)
!
    external :: mmult, mprec, smult
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!
!   allocate memory for the a and s matrices
!
    allocate (a(n,n), s(n,n))
!
!   build a symmetric, positive definite metric matrix 
!
    call random_number(a)
    s = matmul(transpose(a),a)
!
!   build a symmetric matrix:
!
    a = 0.0_dp
    do i = 1, n
      a(i,i) = real(i,kind=dp) + 1.0_dp
      do j = 1, i - 1
        a(j,i) = 1.0_dp/real(i+j,kind=dp)
        a(i,j) = a(j,i)
      end do
    end do
!
!   if required, solve the problem with a dense lapack routine:
!
    if (check_lapack) then
      allocate (a_copy(n,n), s_copy(n,n), w(n))
      a_copy = a
      s_copy = s
      call dsygv(1,'v','l',n,a_copy,n,s_copy,n,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsygv(1,'v','l',n,a_copy,n,s_copy,n,w,work,lwork,info)
!
!     write the results on file for comparison:
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, w(i)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (a_copy(1,i) .lt. 0.0_dp) a_copy(:,i) = - a_copy(:,i)
        write(10,'(10f12.6)') a_copy(:,i)
        write(10,*)
      end do
      close (10)
    end if
!
!   allocate and gather the diagonal:
!
    allocate (diagonal(n))
    do i = 1, n
      diagonal(i) = a(i,i) - s(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvalues and eigenvectors:
!
    allocate (eig(n), evec(n,n_eig))
!
!   compute a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call lobpcg:
!
    call lobpcg_driver(.true.,.true.,n,n_want,n_eig,itmax,tol,0.0_dp,mmult, &
                       mprec,smult,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'lobcpg.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector
!
    call guess_evec(4,n,n_eig,diagonal,evec)
!
!   call the davidson driver:
!
    call gen_david_driver(.true.,n,n_want,n_eig,itmax,tol,m_max,0.0_dp,mmult, &
                         mprec,smult,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'davidson.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_geneig
!
  subroutine test_caslr(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib!, only : caslr_driver, caslr_eff_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   this subroutine builds a generalized eigenvalue problem such as the one encountered
!   in casscf linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical  :: ok
    integer  :: i, j, lwork, info
    integer  :: n2, n_eig
    real(dp) :: sqrttwo, lw(1)
    real(dp), allocatable :: evec(:,:), eig(:), diagonal(:)
    real(dp), allocatable :: work(:), w(:)
!
    external apbvec, ambvec, spdvec, smdvec, lrprec_1, lrprec_2, lambdavec, omegavec, lrprec_2f
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!   the actual size of the problem is 2*n:
!
    n2 = 2 * n
!
!   allocate memory for the a, b, apb, amb, sigma, delta, spd, smd matrices:
!
    allocate (aa(n,n), bb(n,n), apb(n,n), amb(n,n), sigma(n,n), delta(n,n), spd(n,n), smd(n,n))
!
!   build a positive definite, symmetric apb, amb, sigma matrices:
!  
    do i = 1, n
      apb(i,i) = 5.0_dp + real(i,kind=dp)
      do j = 1, i - 1
        apb(j,i) = 1.0_dp/real(i+j,kind=dp)
        apb(i,j) = apb(j,i)
      end do 
    end do
    do i = 1, n
      amb(i,i) = 2.0_dp + dble(i)
      do j = 1, i - 1
        amb(j,i) = 0.2_dp/dble(i+j)
        amb(i,j) = amb(j,i)
      end do 
    end do
!  
    call random_number(sigma)
    delta = matmul(transpose(sigma),sigma)
    sigma = delta
    do i = 1, n
      sigma(i,i) = sigma(i,i) + 100.0_dp
    end do
!  
!   build antisymmetric delta:
!  
    call random_number(delta)
    delta = delta - transpose(delta)
!  
!   build a and b:
!  
    aa = 0.5_dp * (apb + amb)
    bb = 0.5_dp * (apb - amb)
!  
!   build sigma + delta and sigma - delta
!  
    spd = sigma + delta
    smd = sigma - delta 
!
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (a(n2,n2), s(n2,n2), w(n2), lambdafull(n2,n2), omegafull(n2,n2))
      a(1:n,   1:n)    = aa
      a(n+1:n2,n+1:n2) = aa
      a(1:n,   n+1:n2) = bb
      a(n+1:n2,1:n)    = bb
      s(1:n,   1:n)    = sigma
      s(n+1:n2,n+1:n2) = - sigma
      s(1:n,   n+1:n2) = delta
      s(n+1:n2,1:n)    = - delta
!
      lambdafull = a
      omegafull  = s
!
!   if required, solve the generalized eigenvalue problem with a dense
!   lapack routine:
!
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,work,lwork,info)
!
!     write the results on a file for comparison.
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, 1.0_dp/w(n2-i+1)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (s(1,n2-i+1) .lt. 0.0_dp) s(:,n2-i+1) = - s(:,n2-i+1)
        write(10,'(10f12.6)') s(:,n2-i+1)
        write(10,*)
      end do
      close (10)
!
!     free the memory:
!
      deallocate (a, s, w, work)
    end if
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n
      diagonal(i) = aa(i,i) - sigma(i,i)
    enddo
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvectors and eigenvalues:
!
    allocate (evec(n2,n_eig), eig(n_eig))
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(1,n2,n_eig,diagonal,evec)
!
!   call the traditional solver:
!
    write(6,*) ' traditional implementation'
    write(6,*)
    call caslr_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                      spdvec,smdvec,lrprec_1,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    sqrttwo = sqrt(2.0_dp)
    open (unit = 10, file = 'caslr.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/sqrttwo
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the traditional solver using the helmich-paris algorithm:
!
    i_alg = 1
    write(6,*) ' helmich-paris implementation'
    write(6,*)
    call caslr_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                      spdvec,smdvec,lrprec_1,eig,evec,ok)
    i_alg = 0
!
!   write the converged results on file for comparison:
!
    sqrttwo = sqrt(2.0_dp)
    open (unit = 10, file = 'cashp.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/sqrttwo
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the modified solver:
!
    write(6,*) ' new implementation'
    write(6,*)
    call caslr_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                          spdvec,smdvec,lrprec_2,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'caslr_eff.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/2.0_dp
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_caslr
!
  subroutine test_caslr_newcomplex(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    integer              :: pwhat
    real(dp), intent(in) :: tol
!
!   this subroutine builds a generalized eigenvalue problem such as the one encountered
!   in casscf linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical               :: ok
    integer               :: i, j, lwork, info, i_eig
    integer               :: n2, n4, n_eig
    real(dp)              :: sqrttwo, lw(1)
    real(dp), allocatable :: evecre(:,:), evecim(:,:), eig(:), diagonal(:), evec(:,:)
    real(dp), allocatable :: work(:), w(:), wre(:), wim(:)
!
!   debug
!
    real(dp), allocatable :: lambda4(:,:), omega4(:,:) 
    real(dp), allocatable :: prod(:), prod2(:)
!
    external apbvec_newcomplex, ambvec_newcomplex, spdvec_newcomplex, smdvec_newcomplex, &
             lrprec_1_newcomplex, lrprec_2_newcomplex, lambdavec, omegavec, lrprec_2f
!
    1000 format(t3,' eigenvalue # ',i6,': ',f15.6,/,t3,' eigenvector: ')
!
    write(6,*) 'NEW COMPLEX IMPLEMENTATION: MAIN IS ENTERED'
!
!   the actual size of the problem is 2*n:
!
    n2 = 2 * n
    n4 = 4 * n
!
!   allocate memory for the a, b, apb, amb, sigma, delta, spd, smd matrices:
!
    allocate (aare(n,n), aaim(n,n), bbre(n,n), bbim(n,n), &         !Ar, Ai, Br and Bi matrices 
              apbre(n,n), apbim(n,n), ambre(n,n), ambim(n,n), spdre(n,n), spdim(n,n), smdre(n,n), smdim(n,n), & 
              sigmare(n,n), sigmaim(n,n), deltare(n,n), deltaim(n,n)) 
!
!   debug
!
    allocate(lambda4(n4,n4), omega4(n4,n4), prod(n4), prod2(n4))
    prod    = 0.0_dp
    prod2   = 0.0_dp
    lambda4 = 0.0_dp
    omega4  = 0.0_dp
!
!   build a positive definite, hermitian apb, amb matrices:
!   hermitian matrices: Re(sym), Im(antisym)
!  
    apbre = 0.0_dp
    apbim = 0.0_dp
    ambim = 0.0_dp
    ambim = 0.0_dp
    do i = 1, n
      apbre(i,i) = 5.0_dp + real(i,kind=dp)
      do j = 1, i - 1
        apbre(j,i) = 1.0_dp/real(i+j,kind=dp)
        apbre(i,j) = apbre(j,i)
        apbim(j,i) = 1.0_dp/real(i+j,kind=dp)
        apbim(i,j) = -apbim(j,i)
      end do 
    end do
    do i = 1, n
      ambre(i,i) = 2.0_dp + dble(i)
      do j = 1, i - 1
        ambre(j,i) = 0.2_dp/dble(i+j)
        ambre(i,j) = ambre(j,i)
        ambim(j,i) = 0.2_dp/dble(i+j)
        ambim(i,j) = -ambim(j,i)
      end do 
    end do
!
!   build sigma: hermitian
!
    sigmare = 0.0_dp
    sigmaim = 0.0_dp
    deltare = 0.0_dp
    deltaim = 0.0_dp
    call random_number(sigmare)
    deltare = matmul(transpose(sigmare),sigmare)
    sigmare = deltare
    call random_number(sigmaim)
    !deltaim = matmul(transpose(sigmaim),sigmaim)
    sigmaim = (sigmaim - transpose(sigmaim))/2.0_dp
    do i = 1, n
      sigmare(i,i) = sigmare(i,i) + 30.0_dp
    end do
!  
!   build antisymmetric delta:
!   antisymmetric matrices: Re(antisym), Im(sym)
!  
    call random_number(deltare)
    deltare = (deltare - transpose(deltare))/2.0_dp
    call random_number(deltaim)
    deltaim = (deltaim - transpose(deltaim))/2.0_dp
!  
!   build a and b:
!   A is hermitian, B is symmetric
!  
    aare = 0.0_dp
    aaim = 0.0_dp
    bbre = 0.0_dp
    bbim = 0.0_dp
!
    aare = 0.5_dp * (apbre + ambre)
    bbre = 0.5_dp * (apbre - ambre)
    aaim = 0.5_dp * (apbim + ambim)
    bbim = 0.5_dp * (apbim - ambim)
    do i = 1, n
      do j = 1, i
        if (i .ne. j) bbim(i,j) = - bbim(i,j)
      end do
    end do
    apbim = aaim + bbim
    ambim = aaim - bbim
!  
!   build sigma + delta and sigma - delta
!  
    spdre = sigmare + deltare
    smdre = sigmare - deltare 
    spdim = sigmaim + deltaim
    smdim = sigmaim - deltaim 
!
!   to simplify the problem with an SCF-like metric
!
!     sigmare = 0.0_dp
!     sigmaim = 0.0_dp
!     do i = 1, n
!       sigmare(i,i) = 1.0_dp
!       sigmaim(i,i) = 1.0_dp  
!     end do
!     deltare = 1.0_dp 
!     deltaim = 1.0_dp 
!
!   debug: print the matrices
! 
!dc    write(6,*) '(A+B)re: sym'
!dc    call prtmat(6,n,apbre,1)
!dc    write(6,*) '(A+B)im: antisym'
!dc    call prtmat(6,n,apbim,1)
!dc    write(6,*) '(A-B)re: sym'
!dc    call prtmat(6,n,ambre,1)
!dc    write(6,*) '(A-B)im: antisym'
!dc    call prtmat(6,n,ambim,1)
!dc    write(6,*) 'sigmare: sym'
!dc    call prtmat(6,n,sigmare,1)
!dc    write(6,*) 'sigmaim: antisym'
!dc    call prtmat(6,n,sigmaim,1)
!dc    write(6,*) 'deltare: antisym'
!dc    call prtmat(6,n,deltare,1)
!dc    write(6,*) 'deltaim: sym'
!dc    call prtmat(6,n,deltaim,1)
!dc    write(6,*) 'A re: sym'
!dc    call prtmat(6,n,aare,1)
!dc    write(6,*) 'A im: antisym'
!dc    call prtmat(6,n,aaim,1)
!dc    write(6,*) 'B re: sym'
!dc    call prtmat(6,n,bbre,1)
!dc    write(6,*) 'B im: sym'
!dc    call prtmat(6,n,bbim,1)
!
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (are(n2,n2), aim(n2,n2), sre(n2,n2), sim(n2,n2), w(n4), wre(n2), wim(n2), lambdafull(n4,n4), omegafull(n4,n4))
!     real part 
      are(1:n,   1:n)    = aare
      are(n+1:n2,n+1:n2) = aare
      are(1:n,   n+1:n2) = bbre
      are(n+1:n2,1:n)    = bbre
      sre(1:n,   1:n)    = sigmare
      sre(n+1:n2,n+1:n2) = - sigmare
      sre(1:n,   n+1:n2) = deltare
      sre(n+1:n2,1:n)    = - deltare
!     imaginary part       
      aim(1:n,   1:n)    = aaim
      aim(n+1:n2,n+1:n2) = - aaim
      aim(1:n,   n+1:n2) = bbim
      aim(n+1:n2,1:n)    = - bbim
      sim(1:n,   1:n)    = sigmaim
      sim(n+1:n2,n+1:n2) = sigmaim
      sim(1:n,   n+1:n2) = deltaim
      sim(n+1:n2,1:n)    = deltaim
!
!dc      write(6,*) 'are'
!dc      call prtmat(n2,n2,are,1)
!dc      write(6,*) 'aim'
!dc      call prtmat(n2,n2,aim,1)
!dc      write(6,*) 'sigmare: sym'
!dc      call prtmat(n,n,sigmare,1)
!dc      write(6,*) 'sigmaim: antisym'
!dc      call prtmat(n,n,sigmaim,1)
!dc      write(6,*) 'deltare: antisym'
!dc      call prtmat(n,n,deltare,1)
!dc      write(6,*) 'deltaim: sym'
!dc      call prtmat(n,n,deltaim,1)
!dc      write(6,*) 'sre'
!dc      call prtmat(n2,n2,sre,1)
!dc      write(6,*) 'sim'
!dc      call prtmat(n2,n2,sim,1)
!
!     build the 4nx4n matrices lambda and omega 
!
      lambda4(1:n2,1:n2)       = are 
      lambda4(1:n2,n2+1:n4)    = - aim 
      lambda4(n2+1:n4,1:n2)    = aim
      lambda4(n2+1:n4,n2+1:n4) = are 
      omega4(1:n2,1:n2)        = sre 
      omega4(1:n2,n2+1:n4)     = - sim 
      omega4(n2+1:n4,1:n2)     = sim
      omega4(n2+1:n4,n2+1:n4)  = sre 
!
      lambdafull = lambda4
      omegafull  = omega4
!      
!dc      write(6,*) 'Lambda4'
!dc      call prtmat(n4,n4,lambda4,1)
!dc      write(6,*) 'Omega4'
!dc      call prtmat(n4,n4,omega4,1)
!
!     if required, solve the generalized eigenvalue problem with a dense
!     lapack routine:
!
       call dsygv(1,'V','L',n4,omega4,n4,lambda4,n4,w,lw,-1,info)
       print *, 'info1 lapack', info
       lwork = int(lw(1))
       allocate (work(lwork))
       call dsygv(1,'V','L',n4,omega4,n4,lambda4,n4,w,work,lwork,info)
       print *, 'info2 lapack', info
!
!     write the results on a file for comparison.
!
      open (unit = 10, file = 'lapack_newcomplex.txt', form = 'formatted', access = 'sequential')
      do i = 1, n4 !n4!n_want
        write(10,1000) i, 1.0_dp/w(n4-i+1)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (omega4(1,n4-i+1) .lt. 0.0_dp) omega4(:,n4-i+1) = - omega4(:,n4-i+1)
        write(10,'(10f15.6)') omega4(:,n4-i+1)
        write(10,*)
      end do
      close (10)
!
!     free the memory:
!
      deallocate (omega4, lambda4, aim, sim, w, wre, wim, work)
    end if
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n2))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n2
      diagonal(i) = are(i,i) - sre(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
    n_eig = n_want*2 !+ 10
!    n_eig = n4
!
!   allocate memory for the eigenvectors and eigenvalues:
!
    allocate (evecre(n2,n_eig), evecim(n2,n_eig), eig(n_eig), evec(n4,n_eig))
!
!   decide which solver you need: traditional or smogd (both complex) 
!   
    write(6,*)
    write(6,*) 'CHOOSE YOUR FIGHTER:'
    write(6,*)
    write(6,*) '1 - traditional complex solver'
    write(6,*)
    write(6,*) '2 - complex smogd'
    write(6,*)
    read(5,*) pwhat
!
    if (pwhat.eq.1) then
!
!     make a guess for the eigenvector (see guess_evec_newcomplex for more information)
!
!     1: random guess * 0.01_dp
!     2: minloc for both evecre and evecim     
!     3: maxloc for both evecre and evecim
! 
      !call guess_evec_newcomplex(1,n2,n_eig,diagonal,evecre,evecim)
      call guess_evec(4,n4,n_eig,diagonal,evec)
!      do i = 1, n_eig, 2
!        evec(1:n2,i+1)    = -evec(n2+1:n4,i)
!        evec(1:n,i+1)     = -evec(n2+1:n*3,i)
!        evec(n+1:n2,i+1)  =  evec(n*3+1:n4,i)
!        evec(n2+1:n4,i+1) =  evec(1:n2,i)
!      end do
!
!     evec structure: X = (Yr Zr Yi -Zi) 
!
      evecre = evec(1:n2,:)
      evecim = evec(n2+1:n4,:)
!
!     call the traditional fake solver which uses the whole expansion vectors and matrices:
!
!2x2      write(6,*)
!2x2      write(6,*) ' traditional implementation - complex fake Olsen'
!2x2      write(6,*)
!2x2      !call caslr_newcomplex_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec_newcomplex, &
!2x2      !                             ambvec_newcomplex,spdvec_newcomplex,smdvec_newcomplex, &
!2x2      !                             lrprec_1_newcomplex,eig,evecre,evecim,ok)
!2x2      call caslr_2x2_driver(.true.,n2,n4,n_want,n_eig,itmax,tol,m_max, &
!2x2                            lambdavec,omegavec,lrprec_2f,eig,evec,ok)
!2x2      !call caslr_4x4_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec, &
!2x2      !                             ambvec,spdvec,smdvec, &
!2x2      !                             lrprec_1,eig,evec,ok)
!2x2!
!2x2!     write the converged results on file for comparison:
!2x2!
!2x2      sqrttwo = sqrt(2.0_dp)
!2x2      open (unit = 10, file = 'caslr_fake_cOlsen.txt', form = 'formatted', access = 'sequential')
!2x2      do i = 1, n_want
!2x2        write(10,1000) i, eig(i)
!2x2        !if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
!2x2        write(10,'(10f12.6)') evec(:,i)/sqrttwo
!2x2        write(10,*)
!2x2      end do
!2x2      close (10)
!      
!     test A x=lambda S x
!
!test      do i_eig = 1, n_want
!test        write(6,*) 'TEST AX=LAMBDASX', eig(i_eig)
!test        write(6,*) 'eigenval n_', i_eig
!test        call dgemv('n',n4,n4,1.0_dp,lambdafull,n4,evec(:,i_eig),1,0.0_dp,prod,1)
!test        call dgemv('n',n4,n4,1.0_dp,omegafull,n4,evec(:,i_eig),1,0.0_dp,prod2,1)
!test        !write(6,'(2a15)') 'AX', 'LAMBDA SX'
!test        write(6,'(2a15)') 'AX - LAMBDA SX'
!test        do i = 1, n4
!test          !write(6,'(2f15.8)') prod(i), eig(i_eig)*prod2(i)
!test          write(6,'(f15.8)') prod(i) - eig(i_eig)*prod2(i)
!test        end do
!test      end do
!
!     call the traditional solver:
!
      write(6,*)
      write(6,*) ' traditional implementation - complex Olsen'
      write(6,*)
!
      call caslr_newcomplex_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec_newcomplex, &
                                   ambvec_newcomplex,spdvec_newcomplex,smdvec_newcomplex, &
                                   lrprec_1_newcomplex,eig,evecre,evecim,ok)
!
!     write the converged results on file for comparison:
!
      sqrttwo = sqrt(2.0_dp)
      open (unit = 10, file = 'caslr_trad_newcomplex.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, eig(i)
        if (evecre(1,i) .lt. 0.0_dp) evecre(:,i) = - evecre(:,i)
        if (evecim(1,i) .lt. 0.0_dp) evecim(:,i) = - evecim(:,i)
        write(10,'(10f12.6)') evecre(:,i)/sqrttwo
        write(10,'(10f12.6)') evecim(:,i)/sqrttwo
        write(10,*)
      end do
      close (10)
!      
!     test A x=lambda S x
!
!test      evec(1:n2,:)    = evecre
!test      evec(n2+1:n4,:) = evecim 
!test      do i_eig = 1, n_want
!test        write(6,*) 'TEST AX=LAMBDASX', eig(i_eig)
!test        write(6,*) 'eigenval n_', i_eig
!test        call dgemv('n',n4,n4,1.0_dp,lambdafull,n4,evec(:,i_eig),1,0.0_dp,prod,1)
!test        call dgemv('n',n4,n4,1.0_dp,omegafull,n4,evec(:,i_eig),1,0.0_dp,prod2,1)
!test        !write(6,'(2a15)') 'AX', 'LAMBDA SX'
!test        write(6,'(a15)') 'AX - LAMBDA SX'
!test        do i = 1, n4
!test          !write(6,'(2f15.8)') prod(i), eig(i_eig)*prod2(i)
!test          write(6,'(f15.8)') prod(i) - eig(i_eig)*prod2(i)
!test        end do
!test      end do
!
    else
!
!     make a guess for the eigenvector (see guess_evec_newcomplex for more information)
!
      call guess_evec(1,n4,n_eig,diagonal,evec)
!
!     evec structure: X = (Yr Zr Yi -Zi) 
!
      evecre = evec(1:n2,:)
      evecim = evec(n2+1:n4,:)
!      
!     1: minloc for both evecre and evecim     
!
!
!     call the traditional solver:
!
      write(6,*)
      write(6,*) ' traditional implementation - complex Olsen'
      write(6,*)
!
      call caslr_newcomplex_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec_newcomplex, &
                                   ambvec_newcomplex,spdvec_newcomplex,smdvec_newcomplex, &
                                   lrprec_1_newcomplex,eig,evecre,evecim,ok)
!
!     write the converged results on file for comparison:
!
      sqrttwo = sqrt(2.0_dp)
      open (unit = 10, file = 'caslr_trad_newcomplex.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, eig(i)
        if (evecre(1,i) .lt. 0.0_dp) evecre(:,i) = - evecre(:,i)
        if (evecim(1,i) .lt. 0.0_dp) evecim(:,i) = - evecim(:,i)
        write(10,'(10f12.6)') evecre(:,i)/sqrttwo
        write(10,'(10f12.6)') evecim(:,i)/sqrttwo
        write(10,*)
      end do
      close (10)
!
!     make a guess for the eigenvector
!
      call guess_evec(1,n4,n_eig,diagonal,evec)
!
!     evec structure: X = (Yr Zr Yi -Zi) 
!
      evecre = evec(1:n2,:)
      evecim = evec(n2+1:n4,:)
!
!      call guess_evec_newcomplex(1,n2,n_eig,diagonal,evecre,evecim)
!
!     call the modified solver:
!
      write(6,*)
      write(6,*) ' efficient implementation - complex SMO-GD'
      write(6,*)
!
      call caslr_newcomplex_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec_newcomplex,ambvec_newcomplex, &
                                       spdvec_newcomplex,smdvec_newcomplex, & 
                                       lrprec_2_newcomplex,eig,evecre,evecim,ok)
!
!     write the converged results on file for comparison:
!
      sqrttwo = sqrt(2.0_dp)
      open (unit = 10, file = 'caslr_newcomplex_eff.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, eig(i)
        if (evecre(1,i) .lt. 0.0_dp) evecre(:,i) = - evecre(:,i)
        if (evecim(1,i) .lt. 0.0_dp) evecim(:,i) = - evecim(:,i)
        write(10,'(10f12.6)') evecre(:,i)/sqrttwo
        write(10,'(10f12.6)') evecim(:,i)/sqrttwo
        write(10,*)
      end do
      close (10)
    end if
!
    return
  end subroutine test_caslr_newcomplex
!
  subroutine test_caslr_complex(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib!, only : caslr_driver, caslr_complex_eff_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   this subroutine builds a generalized eigenvalue problem such as the one encountered
!   in casscf linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical                   :: ok
    integer                   :: i, j, lwork, info, cwhat, liwork, lrwork
    integer                   :: n2, n_eig
    integer,     allocatable  :: iwork(:)
    complex(dp)               :: sqrttwo
    real(dp)                  :: lw(1)
    real(dp),    allocatable  :: diagonal(:), eig(:)
    real(dp),    dimension(4) :: vecr
    real(dp),    allocatable  :: rwork(:), w(:)
    complex(dp)               :: csqrttwo
    complex(dp), parameter    :: cone=(1.0_dp,0.0_dp), ctwo=(2.0_dp,0.0_dp)
    complex(dp), allocatable  :: cw(:)
    integer,     dimension(4) :: vec
    complex(dp), allocatable  :: evec(:,:)
    complex(dp), allocatable  :: cwork(:)
    complex(dp), allocatable  :: ceig(:)
    real(dp),    allocatable  :: sigma_re(:,:), sigma_im(:,:), delta_re(:,:), delta_im(:,:)
!   to debug: dc
    complex(dp), allocatable  :: abbacopy(:,:), sddscopy(:,:), cprod(:), cprod2(:)
!
    external apbvec_complex, ambvec_complex, spdvec_complex, smdvec_complex, lrprec_1_complex, lrprec_2_complex
!
    1000 format(t3,' eigenvalue # ',i6,': ',f20.10,/,t3,' eigenvector: ')
!
!   the actual size of the problem is 2*n:
!
    n2 = 2 * n
!
!   allocate memory for the a, b, sigma, delta, matrices:
!
    allocate (caa(n,n), cbb(n,n), capb(n,n), camb(n,n), csigma(n,n), cdelta(n,n), abbacopy(n2,n2), sddscopy(n2,n2), &
              cprod(n2), cprod2(n2), sigma_re(n,n), sigma_im(n,n), delta_re(n,n), delta_im(n,n))
!
!   build a positive definite, symmetric A+B, A-B, Sigma, Delta matrices:
!
    sigma_re = 0.0_dp
    sigma_im = 0.0_dp
    delta_re = 0.0_dp
    delta_im = 0.0_dp
    caa      = (0.0_dp,0.0_dp)
    cbb      = (0.0_dp,0.0_dp)
    capb     = (0.0_dp,0.0_dp)
    camb     = (0.0_dp,0.0_dp) 
    csigma   = (0.0_dp,0.0_dp)
    cdelta   = (0.0_dp,0.0_dp)
    abbacopy = (0.0_dp,0.0_dp)
    sddscopy = (0.0_dp,0.0_dp)
!
    do i = 1, n
      capb(i,i) = dcmplx(i+5.0_dp,0.0_dp)
      do j = 1, i - 1
        capb(j,i) = dcmplx(1.0_dp/(i+j),1.0_dp/(i+j))
        capb(i,j) = conjg(capb(j,i))
      end do 
    end do
    do i = 1, n
      camb(i,i) = dcmplx(i+2.0_dp,0.0_dp)
      !write(6,*) camb(i,i)
      do j = 1, i - 1
        camb(j,i) = dcmplx(0.2_dp/(i+j),0.2_dp/(i+j))
        camb(i,j) = conjg(camb(j,i))
      end do 
    end do
    print *, 'capb'
    call prtmat_complex(4,n,capb,3)!    3 is the free-format to correctly print a complex matrix
    print *, 'camb'
    call prtmat_complex(4,n,camb,3)!    3 is the free-format to correctly print a complex matrix
!    capb = (capb + conjg(transpose(capb)))
!    camb = (camb + conjg(transpose(camb)))
!  
!   build A and B, complex, positive definite matrices:
!  
    caa = (capb + camb) * (0.5_dp,0.0_dp)         !dcmplx(0.5_dp,0.5_dp)
    cbb = (capb - camb) * (0.5_dp,0.0_dp)          !dcmplx(0.5_dp,0.5_dp)
!    caa = (caa + conjg(transpose(caa)))/2.0_dp  
!    cbb = (cbb + transpose(cbb))/2.0_dp
    do i = 1, n
      do j = 1, i 
        if (i .ne. j) cbb(i,j) = conjg(cbb(i,j))
      end do
    end do
    print *, 'caa'
    call prtmat_complex(4,n,caa,3)!    3 is the free-format to correctly print a complex matrix
    print *, 'cbb'
    call prtmat_complex(4,n,cbb,3)!    3 is the free-format to correctly print a complex matrix
    !stop
!
!   caa is the complex and hermitian matrix A 
!
!    call random_number(vecr)
!    vec = floor(101*vecr)
!    if (mod(vec(4),2) .ne. 1) vec(4)=vec(4)+1
!    call zlarnv(4,vec,n*n,caa)
!    caa = (caa + conjg(transpose(caa)))/2.d0  
!
!   cbb is the complex and symmetric matrix B 
!
!    call random_number(vecr)
!    vec = floor(101*vecr)
!    if (mod(vec(4),2) .ne. 1) vec(4)=vec(4)+1
!    call zlarnv(4,vec,n*n,cbb)
!    cbb = (cbb + transpose(cbb))/2.d0
!
!  csigma is the complex, hermitian and positive definite matrix sigma
!
!  use the following lines to generate a random complex matrix
!
    !call random_number(vecr)
    !vec = floor(101*vecr,dp)
    !write(6,*) 'vec', vec
    !if (mod(vec(4),2) .ne. 1) vec(4)=vec(4)+1
    !call zlarnv(4,vec,n*n,csigma)
!
!   to obtain sigma complex random 
!
    call random_number(sigma_re)
    call random_number(sigma_im)
    sigma_re = sigma_re - 0.5_dp
    sigma_im = sigma_im - 0.5_dp
   !write(6,*) 'sigmare', sigma_re
   !write(6,*) 'sigmaim', sigma_im 
    csigma = dcmplx(sigma_re,sigma_im)

!    print *, 'csigma', shape(csigma)
!    call prtmat_complex(4,n,csigma,3)!    3 is the free-format to correctly print a complex matrix
    cdelta = matmul(conjg(transpose(csigma)),csigma)
    csigma = cdelta
    do i = 1, n
      csigma(i,i) = csigma(i,i) + (20.0_dp,0.0_dp)
    end do  
!
!   new sigma&delta: simpler but complex
!
!    do i = 1, n
!      csigma(i,i) = dcmplx(i+3.0_dp,0.0_dp)
!      do j = 1, i - 1
!        csigma(j,i) = dcmplx(1.0_dp/(i+j),0.0_dp/(i+j))
!        csigma(i,j) = conjg(csigma(j,i))
!      end do 
!    end do
!    do i = 1, n
!      cdelta(i,i) = (0.0_dp,0.0_dp)
!      do j = 1, i - 1
!        cdelta(j,i) = (1.0_dp, 1.0_dp)
!        cdelta(i,j) = (-1.0_dp, -1.0_dp)
!      end do
!    end do
!
!  cdelta is the complex, antisymmetric matrix delta
!  use the following lines to generate a random complex matrix
!
!    call random_number(vecr)
!    vec = floor(101*vecr,dp)
!    if (mod(vec(4),2) .ne. 1) vec(4)=vec(4)+1
!    call zlarnv(4,vec,n*n,cdelta)
!
    cdelta = (0.0_dp,0.0_dp)
    call random_number(delta_re)
    call random_number(delta_im)
    delta_re = delta_re - 0.5_dp
    delta_im = delta_im - 0.5_dp
    cdelta = dcmplx(delta_re,delta_im)
    cdelta = (cdelta - transpose(cdelta))/(2.0_dp,0.0_dp) !dcmplx(2.0_dp,2.0_dp)  
!    stop
!
!   to simplify the problem with an SCF-like metric
!
     csigma = (0.0_dp,0.0_dp)
     do i = 1, n
       csigma(i,i) = (1.0_dp,0.0_dp)
     end do
     cdelta = (0.0_dp,0.0_dp)
!
    !print *, 'csigma'
    !call prtmat_complex(4,n,csigma,3)!    3 is the free-format to correctly print a complex matrix
    !print *, 'cdelta'
    !call prtmat_complex(4,n,cdelta,3)!    3 is the free-format to correctly print a complex matrix
!
!    To obtain a real problem:    
!
!      caa    = (caa + conjg(caa))/2.0_dp 
!      cbb    = (cbb + conjg(cbb))/2.0_dp 
!      csigma = (csigma + conjg(csigma))/2.0_dp 
!      cdelta = (cdelta + conjg(cdelta))/2.0_dp 
!
!    write(6,*) 'csigma', csigma
!    write(6,*) 'cdelta', cdelta
!
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (ca(n2,n2), cs(n2,n2), cw(n2), w(n2))
      ca(1:n,   1:n)    = caa
      ca(n+1:n2,n+1:n2) = conjg(caa)
      ca(1:n,   n+1:n2) = cbb
      ca(n+1:n2,1:n)    = conjg(cbb)
      cs(1:n,   1:n)    = csigma
      cs(n+1:n2,n+1:n2) = - conjg(csigma)
      cs(1:n,   n+1:n2) = cdelta
      cs(n+1:n2,1:n)    = - conjg(cdelta)
!      
!     to debug: whole ca and cs matrices
!
      abbacopy = ca
!      write(6,*) 'ca', ca
      sddscopy = cs
!      write(6,*) 'cs', cs
!
!   To diagonalize csigma
!
    !allocate (cwork(1), rwork(1), iwork(1))
    !!call zheevd('v','l',n2,cs,n2,w,cwork,-1,rwork,-1,iwork,-1,info)
    !call zheevd('v','l',n,csigma,n,w,cwork,-1,rwork,-1,iwork,-1,info)
    !write(6,*) 'info1', info
    !lwork = int(cwork(1))
    !lrwork = int(rwork(1))
    !liwork = int(iwork(1))
    !deallocate(cwork, rwork, iwork)
    !allocate (cwork(lwork), rwork(lrwork), iwork(liwork))
    !!call zheevd('v','l',n2,cs,n2,w,cwork,lwork,rwork,lrwork,iwork,liwork,info)
    !call zheevd('v','l',n,csigma,n,w,cwork,lwork,rwork,lrwork,iwork,liwork,info)
    !deallocate(rwork,iwork)
    !write(6,*) 'info2', info
    !write(6,*) 'Eigenvalues:'
    !write(6,'(f15.8)') w
    !stop
!   

!      cs = (0.0_dp,0.0_dp)
!      do i = 1, n2
!        cs(i,i) = (1.0_dp,0.0_dp)
!      end do
!      print *, 'sigma'
!      call prtmat_complex(4,n2,cs,3)!    3 is the free-format to correctly print a complex matrix
!      print *, '-delta*'
!      call prtmat_complex(4,n2,cs(n+1,1),3)!    3 is the free-format to correctly print a complex matrix
!      print *, 'delta'
!      call prtmat_complex(4,n2,cs(1,n+1),3)!    3 is the free-format to correctly print a complex matrix
!      print *, '-sigma*'
!      call prtmat_complex(4,n2,cs(n+1,n+1),3)!    3 is the free-format to correctly print a complex matrix
!      stop
!
!   if required, solve the generalized eigenvalue problem with a dense
!   lapack routine:
!
      allocate (cwork(1), rwork(1))
      call zhegv(1,'V','L',n2,cs,n2,ca,n2,w,cwork,-1,rwork,info)
      write(6,*) 'info1 densa', info
      lwork = int(cwork(1))
      deallocate (cwork,rwork)
      allocate (cwork(lwork), rwork(3*n2-2))
      call zhegv(1,'V','L',n2,cs,n2,ca,n2,w,cwork,lwork,rwork,info)
      write(6,*) 'info2 densa', info
      !w = real(cw)
!
!     write the results on a file for comparison.
!
      open (unit = 10, file = 'lapack_complex_caslr.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, 1.0_dp/w(n2-i+1)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (abs(cs(1,n2-i+1)) .lt. 0.0_dp) cs(:,n2-i+1) = - cs(:,n2-i+1)
        write(10,'(10f12.6)') cs(:,n2-i+1)
        write(10,*)
      end do
      close (10)
!
!     free the memory:
!
      deallocate (ca, cs, w, cwork, rwork)
!      deallocate (w, cwork, rwork)
    end if
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n
      diagonal(i) = caa(i,i) - csigma(i,i)
    end do
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
    n_eig = n_want 
    !n_eig = n
!
!   allocate memory for the eigenvectors and eigenvalues:
!
    allocate (evec(n2,n_eig), eig(n_eig), ceig(n_eig))
!
!   decide which solver you need: traditional or smogd (both complex) 
!   
    write(6,*)
    write(6,*) 'CHOOSE YOUR FIGHTER:'
    write(6,*)
    write(6,*) '1 - traditional complex solver'
    write(6,*)
    write(6,*) '2 - complex smogd'
    write(6,*)
    read(5,*) cwhat
!
    if (cwhat.eq.1) then
!   
!     make a guess for the eigenvector (see guess_evec for more information)
!  
!     5: minloc, 6: random
      call guess_evec_complex(6,n2,n_eig,diagonal,evec)
!    
!     X = (Y Z*)
!    
      do i = 1, n_eig
        !evec(n+1:n2,i) = conjg(evec(n+1:n2,i))
        !evec(n+1:n2,i) = conjg(evec(1:n,i))
      end do
!  
!     call the traditional solver:
!  
      write(6,*) ' traditional implementation - complex linear response'
      write(6,*)
      call caslr_complex_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec_complex,ambvec_complex, &
                                spdvec_complex,smdvec_complex,lrprec_1_complex,eig,ceig,evec,abbacopy,sddscopy,ok)
!  
!     write the converged results on file for comparison:
!  
      csqrttwo = sqrt(ctwo)
      open (unit = 10, file = 'caslr_trad_complex.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, eig(i)
        if (abs(evec(1,i)) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
        write(10,'(10f12.6)') evec(:,i)/csqrttwo
        write(10,*)
      end do
      close (10)
!
!     test A x=lambda S x
!
      write(6,*) 'TEST AX=LAMBDASX', ceig(1), eig(1)
      call zgemv('n',n2,n2,(1.0_dp,0.0_dp),abbacopy,n2,evec(:,1),1,(0.0_dp,0.0_dp),cprod,1)
      call zgemv('n',n2,n2,(1.0_dp,0.0_dp),sddscopy,n2,evec(:,1),1,(0.0_dp,0.0_dp),cprod2,1)
      write(6,'(4a15)') 'Re(AX)', 'Re(LAMBDA SX)','Im(AX)', 'Im(LAMBDA SX)'
      do i = 1, n2
        write(6,'(4f15.8)') real(cprod(i)), real(ceig(1)*cprod2(i)), aimag(cprod(i)), aimag(ceig(1)*cprod2(i))
      end do
    else 
!
!     make a guess for the eigenvector (see guess_evec for more information)
!    
!     5: minloc, 6: random
      call guess_evec_complex(5,n2,n_eig,diagonal,evec)
!    
!     X = (Y Z*)
!    
      do i = 1, n_eig
        evec(n+1:n2,i) = conjg(evec(n+1:n2,i))
      end do
!    
!     call the modified solver:
!    
      write(6,*) ' new implementation - complex linear response'
      write(6,*)
      call caslr_complex_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec_complex,ambvec_complex, &
                                    spdvec_complex,smdvec_complex,lrprec_2_complex,eig,ceig,evec,ok)
!    
!     write the converged results on file for comparison:
!    
      open (unit = 10, file = 'caslr_complex_eff.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, eig(i)
        if (abs(evec(1,i)) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
        write(10,'(10f12.6)') evec(:,i)/ctwo  !dcmplx(2.0_dp,2.0_dp)
        write(10,*)
      end do
      close (10)
    end if
!
    return
  end subroutine test_caslr_complex
!
  subroutine test_scflr(check_lapack,n,n_want,tol,itmax,m_max)
    use real_precision
    use utils
    use diaglib, only : caslr_driver, caslr_eff_driver
    implicit none
    logical,  intent(in) :: check_lapack
    integer,  intent(in) :: n, n_want, itmax, m_max
    real(dp), intent(in) :: tol
!
!   this subroutine builds a generalized eigenvalue problem such as the one encountered
!   in hf/dft linear response theory, and then uses lapack and iterative routines to
!   solve it.
!
    logical  :: ok
    integer  :: i, j, lwork, info
    integer  :: n2, n_eig
    real(dp) :: sqrttwo, lw(1)
    real(dp), allocatable :: evec(:,:), eig(:), diagonal(:)
    real(dp), allocatable :: work(:), w(:)
!
    external apbvec, ambvec, spdvec, smdvec, lrprec_1, lrprec_2
!
    1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
!   the actual size of the problem is 2*n:
!
    n2 = 2 * n
!
!   set the flag that avoids multiplications with the metric
!
    tdscf = .true.
!
!   allocate memory for the a, b, apb, amb, sigma, delta, spd, smd matrices:
!
    allocate (aa(n,n), bb(n,n), apb(n,n), amb(n,n), sigma(n,n), delta(n,n), spd(n,n), smd(n,n))
!
!   build a positive definite, symmetric apb, amb matrices:
!  
    do i = 1, n
      apb(i,i) = 5.0_dp + real(i,kind=dp)
      do j = 1, i - 1
        apb(j,i) = 1.0_dp/real(i+j,kind=dp)
        apb(i,j) = apb(j,i)
      end do 
    end do
    do i = 1, n
      amb(i,i) =2.0_dp + dble(i)
      do j = 1, i - 1
        apb(j,i) = 0.2_dp/dble(i+j)
        apb(i,j) = apb(j,i)
      end do 
    end do
!  
!   set sigma to the identity:
!
    sigma = 0.0_dp
    do i = 1, n
      sigma(i,i) = 1.0_dp
    end do
!  
!   set delta to zero:
!  
    delta = 0.0_dp
!  
!   build a and b:
!  
    aa = 0.5_dp * (apb + amb)
    bb = 0.5_dp * (apb - amb)
!  
!   build sigma + delta and sigma - delta
!  
    spd = sigma
    smd = sigma
!
    if (check_lapack) then
!
!     build the complete matrices:
!    
      allocate (a(n2,n2), s(n2,n2), w(n2))
      a(1:n,   1:n)    = aa
      a(n+1:n2,n+1:n2) = aa
      a(1:n,   n+1:n2) = bb
      a(n+1:n2,1:n)    = bb
      s(1:n,   1:n)    = sigma
      s(n+1:n2,n+1:n2) = - sigma
      s(1:n,   n+1:n2) = delta
      s(n+1:n2,1:n)    = - delta
!
!   if required, solve the generalized eigenvalue problem with a dense
!   lapack routine:
!
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork))
      call dsygv(1,'V','L',n2,s,n2,a,n2,w,work,lwork,info)
!
!     write the results on a file for comparison.
!
      open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
      do i = 1, n_want
        write(10,1000) i, 1.0_dp/w(n2-i+1)
!
!       fix the phase so that the first element of the eigenvector is positive.
!
        if (s(1,n2-i+1) .lt. 0.0_dp) s(:,n2-i+1) = - s(:,n2-i+1)
        write(10,'(10f12.6)') s(:,n2-i+1)
        write(10,*)
      end do
      close (10)
!
!     free the memory:
!
      deallocate (a, s, w, work)
    end if
!
!   allocate space for the diagonal:
!
    allocate (diagonal(n))
!
!   gather diag(a) - diag(sigma):
!
    do i = 1, n
      diagonal(i) = aa(i,i) - sigma(i,i)
    enddo
!
!   for better convergence, we seek more eigenpairs and stop the iterations when the
!   required ones are converged.
!
    n_eig = min(2*n_want, n_want + 5)
!
!   allocate memory for the eigenvectors and eigenvalues:
!
    allocate (evec(n2,n_eig), eig(n_eig))
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the traditional solver:
!
    call caslr_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                      spdvec,smdvec,lrprec_1,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    sqrttwo = sqrt(2.0_dp)
    open (unit = 10, file = 'caslr.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/sqrttwo
      write(10,*)
    end do
    close (10)
!
!   make a guess for the eigenvector (see guess_evec for more information)
!
    call guess_evec(4,n2,n_eig,diagonal,evec)
!
!   call the modified solver:
!
    call caslr_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec, &
                          spdvec,smdvec,lrprec_2,eig,evec,ok)
!
!   write the converged results on file for comparison:
!
    open (unit = 10, file = 'caslr_eff.txt', form = 'formatted', access = 'sequential')
    do i = 1, n_want
      write(10,1000) i, eig(i)
      if (evec(1,i) .lt. 0.0_dp) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/2.0_dp
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_scflr
!
  subroutine guess_evec(iwhat,n,m,diagonal,evec)
    use real_precision
    implicit none
    integer,                   intent(in)    :: iwhat, n, m
    real(dp), dimension(n/2),  intent(in)    :: diagonal
    real(dp), dimension(n,m),  intent(inout) :: evec
!
!   guess the eigenvector.
!
    integer               :: i, ipos, n_seed
    integer,  allocatable :: iseed(:)
    real(dp), allocatable :: temp_r(:), temp_i(:)
    logical,  allocatable :: mask(:)
!
!   initialize a random number generator in a predictible way.
!
    call random_seed(size=n_seed)
    allocate (iseed(n_seed))
    iseed = 1
    call random_seed(put=iseed)
    deallocate (iseed)
!
    evec   = 0.0_dp
!
    allocate (mask(n/2))
!
    if (iwhat.eq.1) then
!
!     get the minimum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = 1.0_dp
      end do
    else if (iwhat.eq.2) then
!
!     get the maximum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = maxloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = 1.0_dp
      enddo
    else if (iwhat.eq.3) then
!
!     random vector between 0 and 1
!
      call random_number(evec)
    else if (iwhat.eq.4) then
!
!     random vector between -0.5 and 0.5
!
      call random_number(evec)
      evec = evec - 0.50_dp
    else if (iwhat.eq.5) then
!
      call random_number(evec)
      evec = evec * 0.01_dp
!
!     get the maximum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = maxloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = evec(ipos,i) + 1.0_dp
      enddo
!
    else if (iwhat.eq.6) then
!
      call random_number(evec)
      evec = evec * 0.01_dp
!
!     get the minimum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = evec(ipos,i) + 1.0_dp
      enddo
!
    else if (iwhat.eq.7) then
!
!     get the minimum element of the diagonal
!
      mask = .true.
      allocate(temp_r(n/2), temp_i(n/2))
      do i = 1, m, 2
        !ipos = minloc(diagonal,dim=1,mask=mask)
        !mask(ipos) = .false.   
        !evec(ipos,i)       = 1.0_dp
        !!evec(ipos+n/2,i)   = 1.0_dp
        !!evec(ipos,i+1)     = 1.0_dp
        !evec(ipos+n/2,i+1) = -1.0_dp
        call random_number(temp_r)
        call random_number(temp_i)
        evec(1:n/2,i)     = temp_r
        evec(n/2+1:n,i)   = temp_i
        evec(1:n/2,i+1)   = temp_i
        evec(n/2+1:n,i+1) = -temp_r
      enddo
      deallocate(temp_r, temp_i)
 
    end if
    return
  end subroutine guess_evec
!
  subroutine guess_evec_newcomplex(iwhat,n,m,diagonal,evecre,evecim) 
    use real_precision
    implicit none 
    integer,                  intent(in)    :: iwhat, n, m
    real(dp), dimension(n/2), intent(in)    :: diagonal
    real(dp), dimension(n,m), intent(inout) :: evecre, evecim
!
!   guess the eigenvector.
!
    integer              :: i, ipos, n_seed
    integer, allocatable :: iseed(:)
    logical, allocatable :: mask(:)
!
!   initialize a random number generator in a predictible way.
!
    call random_seed(size=n_seed)
    allocate (iseed(n_seed))
    iseed = 1
    call random_seed(put=iseed)
    deallocate (iseed)
!
    !evec   = 0.0_dp
    evecre = 0.0_dp
    evecim = 0.0_dp
!
    allocate (mask(n/2))
!
    if (iwhat.eq.1) then
!
      call random_number(evecre)
      evecre = evecre * 0.01_dp
!      evecre = evecre - 0.50_dp
      call random_number(evecim)
      evecim = evecim * 0.01_dp
!      evecim = evecim - 0.50_dp
!     
    else if (iwhat.eq.2) then
!     get the minimum element of the diagonal
!     
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evecre(ipos,i) = evecre(ipos,i) + 1.0_dp
        evecim(ipos,i) = evecim(ipos,i) + 1.0_dp
      end do
!
    else if (iwhat.eq.3) then
      mask = .true.
      do i = 1, m
        ipos = maxloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evecre(ipos,i) = evecre(ipos,i) + 1.0_dp
        evecim(ipos,i) = evecim(ipos,i) + 1.0_dp
      end do
    end if 
!
  end subroutine guess_evec_newcomplex
!
  subroutine guess_evec_complex(iwhat,n,m,diagonal,evec)
    use real_precision
    implicit none
    integer,                       intent(in)    :: iwhat, n, m
    real(dp),    dimension(n/2),   intent(in)    :: diagonal
    complex(dp), dimension(n,m),   intent(inout) :: evec
!
!   guess the eigenvector.
!
    integer                 :: i, ipos, n_seed
    integer,  allocatable   :: iseed(:)
    logical,  allocatable   :: mask(:)
    integer,  dimension(4)  :: vec
    real(dp), dimension(4)  :: vecr
!
!   initialize a random number generator in a predictible way.
!
!    call random_seed(size=n_seed)
!    allocate (iseed(n_seed))
!    iseed = 1
!    call random_seed(put=iseed)
!    deallocate (iseed)
!
    evec = (0.0_dp,0.0_dp)
    vecr = 0.0_dp
!
    allocate (mask(n/2))
!
    if (iwhat.eq.5) then
!
      write(6,*) 
      write(6,*) 'GUESS MINLOC'
      write(6,*) 
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.
!        write(6,*) 'evec prima', evec(ipos,i)   
        evec(ipos,i) = evec(ipos,i) + (1.0_dp,0.0_dp)
!        write(6,*) 'evec dopo', evec(ipos,i)   
      end do
      !evec = evec - 0.50_dp
      !evec = evec * 0.1_dp !dcmplx(0.1_dp,0.1_dp)
    else if (iwhat.eq.6) then
      write(6,*) 
      write(6,*) 'GUESS RANDOM'
      write(6,*) 
      call random_number(vecr)
      vec = floor(101*vecr,dp) 
      if (mod(vec(4),2) .ne. 1) vec(4)=vec(4)+1
      call zlarnv(4,vec,n*m,evec)  
      !evec = evec + conjg(evec)
!      evec = evec - (0.50_dp,0.0_dp)
!      evec = (evec + conjg(evec))/2.0_dp
!      
!     qui forse devo scrivere evec(n+1:n2,1:n2...)=conjg(evec(")
!
!     evec = evec * 100_dp
    end if
    return
  end subroutine guess_evec_complex










