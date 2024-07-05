program main
  use real_precision
  use utils
  use diaglib
!
! simple minded test program to call the lobpcg routine
!
  implicit none
  integer  :: n, itmax, m_max, n_want, iwhat
  real(dp) :: tol
!
! initialize:
!
  n      = 2000
  n_want = 10
  tol    = 1.0e-8_dp
  itmax  = 100
  m_max  = 20
  nmult  = 0
  tdscf  = .false.
  i_alg  = 0
!
  call print_header
! 
  write(6,*)
  1000 format(t3,' simple-minded test driver. select ',/, &
              t3,'   1 for symmetric eigenvalue problems ',/, &
              t3,'   2 for symmetric generalized eigenvalue problems ',/, &
              t3,'   3 for linear-response equations (SCF-like)',/, &
              t3,'   4 for linear-response equations (CASSCF-like)',/, &
              t3,'   5 for unsymmetric eigenvalue problems.')
  write(6,1000)
  !read(5,*) iwhat
  write(6,*)
!
 iwhat=5
  if (iwhat.eq.1) then 
    call test_symm(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.2) then 
    call test_geneig(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.3) then 
    call test_scflr(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.4) then 
    call test_caslr(.true.,n,n_want,tol,itmax,m_max)
  else if (iwhat.eq.5) then
    !call test_nonsym(.false.,n,n_want,tol,itmax,m_max,4,.true.,.false.)
    call eval_nonsym()
  else
    write(6,*) ' invalid selection. aborting ...'
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
  subroutine mmult_l(n,m,x,ax)
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
    real(dp), dimension(n,n)                 :: a_t 
    integer :: i, j, icol
!
!   transpose the matrix
!
    a_t = transpose(a)
    nmult = nmult + m
    do icol = 1, m
      ax(:,icol) = matmul(a_t,x(:,icol))
    end do
    return
  end subroutine mmult_l
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
        if (abs(a(i,i)+fac).gt.tol) then
          px(i,icol) = x(i,icol)/(a(i,i) + fac)
        else
          px(i,icol) = x(i,icol)
        end if
      end do
    end do
    return
  end subroutine mprec
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
    yp = 0.0d0
    ym = 0.0d0
!
    do icol = 1, m
      do i = 1, n
        yp(i,icol) = - 1.0d0/(aa(i,i)**2 - fac**2 * sigma(i,i)**2) * (aa(i,i) * xp(i,icol) + fac * sigma(i,i) * xm(i,icol))
        ym(i,icol) = - 1.0d0/(aa(i,i)**2 - fac**2 * sigma(i,i)**2) * (aa(i,i) * xm(i,icol) + fac * sigma(i,i) * xp(i,icol))
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
    yp = 0.0d0
    ym = 0.0d0
!
    do icol = 1, m
      do i = 1, n
        denom = fac * fac * aa(i,i) * aa(i,i) - sigma(i,i) * sigma(i,i)
        denom = 1.0d0/denom
        yp(i,icol) = denom * (fac * aa(i,i) * xp(i,icol) + sigma(i,i) * xm(i,icol))
        ym(i,icol) = denom * (fac * aa(i,i) * xm(i,icol) + sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_2
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_symm
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
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
    use diaglib, only : caslr_driver, caslr_eff_driver
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
    external apbvec, ambvec, spdvec, smdvec, lrprec_1, lrprec_2
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
      amb(i,i) =2.0_dp + dble(i)
      do j = 1, i - 1
        apb(j,i) = 0.2_dp/dble(i+j)
        apb(i,j) = apb(j,i)
      end do 
    end do
!  
    call random_number(sigma)
    delta = matmul(transpose(sigma),sigma)
    sigma = delta
    do i = 1, n
      sigma(i,i) = sigma(i,i) + 1.0_dp
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/2.0_dp
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_caslr
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
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
      if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
      write(10,'(10f12.6)') evec(:,i)/2.0_dp
      write(10,*)
    end do
    close (10)
!
    return
  end subroutine test_scflr
!
  subroutine eval_nonsym()
    use real_precision
    use utils
    implicit none
    integer  :: n, itmax, m_max, n_want, use_mat, iseed, ios, pos, i
    real(dp) :: tol
    logical  :: both, allsvd
    character(len=30) :: line, value_str
!
!   read input file 
!
    open(unit = 10, file = 'diaglib.in', status = 'old', action = 'read', iostat = ios)
    if (ios /= 0) then
        print *, 'Error opening file: diaglib.in does not exist.'
        stop
    endif
!
    do i = 1, 9
!
!     read the line from the file
!
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) then
          print *, 'Error reading file'
          stop
      endif
!

      ! Find the position of the equal sign
      pos = index(line, '=')

      ! Extract the value from the line starting after the equal sign
      value_str = adjustl(line(pos+1:))

      ! Convert the string to an integer
      if (i.eq.1) read(value_str, *) n
      if (i.eq.2) read(value_str, *) n_want
      if (i.eq.3) read(value_str, *) tol
      if (i.eq.4) read(value_str, *) itmax
      if (i.eq.5) read(value_str, *) m_max
      if (i.eq.6) read(value_str, *) use_mat
      if (i.eq.7) read(value_str, *) both
      if (i.eq.8) read(value_str, *) allsvd
      if (i.eq.9) read(value_str, *) iseed
    end do

    close(10)
!
    call test_nonsym(.false.,n,n_want,tol,itmax,m_max,iseed,use_mat,both,allsvd)
  end subroutine
!
  subroutine test_nonsym(check_lapack,n,n_want,tol,itmax,m_max,iseed,use_mat,both,allsvd)
    use real_precision
    use utils
    use diaglib, only : nonsym_driver
    implicit none
    integer, intent(in) :: n, n_want, itmax, m_max, iseed, use_mat
    logical, intent(in) :: check_lapack, both, allsvd
    real(dp), intent(in):: tol
!
!   test matrix
!
    integer                     :: i, j, info, ipiv(n), lwork, seed_size, n_eig
    integer, allocatable        :: seed(:)
    real(dp)                    :: lw(1), zero, one, low, up, fac, dnrm2
    real(dp), allocatable       :: work(:), a_copy(:,:), r(:,:), l(:,:), wr(:), wi(:), diag(:,:), t(:,:), &
                                    p(:,:),eig(:), evec_r(:,:), evec_l(:,:), diagonal(:), expt(:,:), expmt(:,:)
    logical                     :: ok
!
    external :: mmult, mmult_l, mprec
!
    low       = 0.0d0
    up        = 1.d-4
    zero      = 0.d0
    one       = 1.d0
!
!   allocate memory to get the matrix
!
    allocate (a(n,n), diag(n,n), t(n,n), p(n,n), expt(n,n), expmt(n,n), a_copy(n,n))
!    
    diag = zero
!
!   choice of the test matrix (1 = nonsym, 2 = nonsym (pertubation on sym), 3 = sym)
!
    if (use_mat .eq. 1) then
!
!     build a nonsymmetric matrix 
!
!     generate diagonal matrix:   diag(i,i) = 2 + i
!
      forall(i = 1:n) diag(i, i) = real(2 + i, kind=dp)
!
!     generate random matrix:     t = random
!
      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      seed = iseed
      call random_seed(put=seed)
      call random_number(t)
      do i = 1, n
        t(i,i) = t(i,i) + dble(100+i)
      end do
      deallocate(seed)
!
!     ensure positive definite:   p = t^T * t
!
      call dgemm('t','n',n,n,n,one,t,n,t,n,zero,p,n)
      a = p
!
!     get it's invert by using lu decomposition
!
      call dgetrf(n,n,a,n,ipiv,info)
      call dgetri(n,a,n,ipiv,lw,-1,info)
      lwork = int(lw(1))
      allocate(work(lwork))
      call dgetri(n,a,n,ipiv,work,lwork,info)
      deallocate(work)   
!      
!     get final matrix:           a = p * diag * p^-1
!
      call dgemm('n','n',n,n,n,one,p,n,diag,n,zero,t,n)
      call dgemm('n','n',n,n,n,one,t,n,a,n,zero,p,n)
      a = p
!
    else if (use_mat .eq. 2) then
!
!   get a symmetric matrix
!
      do i = 1, n
        do j = i, n
          if (j == i) then
            a(i,j) = i + 1.0d0
          else if (j > i) then
            a(i,j) = 1.0d0 / (dble(i + j))
            a(j,i) = a(i,j)
          end if
        end do
      end do
!
!     generate random matrix :     t = random
!
      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      seed = iseed
      call random_seed(put=seed)
      call random_number(t)
      deallocate(seed)
!
!     adapt to range low:up and set diagonal to zero
!
      t = low + (up-low) * t
!
      do i = 1, n
        t(i,i) = zero
      end do
!
!     introduce perturbation on matrix
!
      a = a + t
!
    else if (use_mat .eq. 3) then 
!
!   get a symmetric matrix
!
      do i = 1, n
        a(i,i) = real(i,kind=dp) + 1.0d0
        do j = i, n
          if (j.ne.i) then
            a(i,j) = 1.0d0 / (dble(i + j))
            a(j,i) = a(i,j)
          end if
        end do
      end do
!
    else if (use_mat .eq. 4) then 
!
!   get an unsymmetric matrix by using matrix exponentials to create a similarity transform.
!
      a = 0.0_dp
      do i = 1, n
        a(i,i) = real(i,kind=dp) + 1.0_dp
      end do
! 
      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      seed = iseed
      call random_seed(put=seed)
      deallocate(seed)
!
      call random_number(t)
!
!     make sure that the matrix has a small norm:
!
      fac  = dnrm2(n*n,t,1)
      fac  = 0.01_dp / fac
      t    = fac * t
!
!     uncomment the following line to obtain a unitary transformation!
!
!     tmat = tmat - transpose(tmat)
!
!     compute exp(tmat):
!
      call matexp(n,t,expt)
!
!     compute exp(- tmat):
!
      t   = - t
      call matexp(n,t,expmt)
!
!     a = exp(-t) a exp(t)
!
      a_copy = matmul(a,expt)
      a      = matmul(expmt,a_copy) 
!
    else if (use_mat .eq. 5) then
!
!     get a symmetric matrix
!
      do i = 1, n
        a(i,i) = real(i,kind=dp) + 1.0d0
        do j = i, n
          if (j.ne.i) then
            a(i,j) = 1.0d0 / (dble(i + j))
            a(j,i) = a(i,j)
          end if
        end do
      end do
!
!     apply similarity transform
! 
      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      seed = iseed
      call random_seed(put=seed)
      deallocate(seed)
!
      call random_number(t)
!
!     make sure that the matrix has a small norm:
!
      fac  = dnrm2(n*n,t,1)
      fac  = 0.01_dp / fac
      t    = fac * t
!
!     uncomment the following line to obtain a unitary transformation!
!
!     tmat = tmat - transpose(tmat)
!
!     compute exp(tmat):
!
      call matexp(n,t,expt)
!
!     compute exp(- tmat):
!
      t   = - t
      call matexp(n,t,expmt)
!
!     a = exp(-t) a exp(t)
!
      a_copy = matmul(a,expt)
      a      = matmul(expmt,a_copy) 
!
    else
      print *, "no valid matrix choice in test_nonsym."
      stop
    end if
!
    deallocate (p, t, diag, expmt, expt, a_copy)
!
!   for better convergence, we seek more eigenpairs and stop the iterations when 
!   the required ones are converged
!
    n_eig =  n_want 
!
!   print some information
!
       
    1000 format(t5,55("-"),/,t3,'   nonsymmetric davidson eigensolver test run',/,t5,55("-"),/)
    1100 format(t5,55("="),/,t3,'   input information',/,t5,55("-"),/, &
                t3,'  dimension of the full space            :   ',i8,/, &
                t3,'  sought number of eigenpairs            :   ',i8,/, &
                t3,'  number of vectors added each iteration :   ',i8,/, &
                t3,'  convergency tolerance of residual norm :   ',d8.2,/, &
                t3,'  maximum iterations                     :   ',i8,/, &
                t3,'  size of expansion space                :   ',i8,/, &
                t3,'  used matrix                            :   ',i8,/, &
                t3,'  seed for matrix generation             :   ',i8,/,&
                t3,'  calculation of left and right pairs    :   ',l8,/,&
                t3,'  use loop over svd in biortho_vs_x      :   ',l8,/,&
                t5,55("="))
    write(6,1000)
    write(6,1100) n,n_want,n_eig,tol,itmax,m_max,use_mat,iseed,both,allsvd
!
!  if required, solve the problem with a dense lapack routine:
!
    if (check_lapack) then
      allocate (a_copy(n,n), wr(n), wi(n), r(n,n), l(n,n))
      a_copy = a
      call dgeev('v','v',n,a_copy,n,wr,wi,l,n,r,n,lw,-1,info)
      lwork = int(lw(1))
      allocate (work(lwork)) 
      call dgeev('v','v',n,a_copy,n,wr,wi,l,n,r,n,work,lwork,info)
      deallocate (work, a_copy, wr, wi, r, l)
!
      if (info.ne.0) then
        print *, "diagonalization of full space failed."
        stop
      end if
!
      call sort_eigenpairs(wr,wi,r,l,n,n,n,n,.true.,1.d-16)
!
      print *
      1200 format(t5,55("-"),/,t3,'   eigenvalues of lapack full space diagonalization',/,t5,55("-"))
      write(6,1200)
      call printMatrix(n_want,1,wr,n)
      print *
      print *
    end if 
!
!   allocate and gather the diagonal
!
    allocate (diagonal(n))
    do i = 1,n
      diagonal(i) = a(i,i)
    end do
    !n_eig =  min(3*n_want, n_want + 5)
!
      1300 format(t5,55("-"),/,t3,'   nonsymmetric davidson results',/,t5,55("-"),/)
      write(6,1300)
!
!   print output file
!
    open (unit = 10, file = 'diaglib.out', form = 'formatted', access = 'sequential')
    write(10,1000)
    write(10,1100) n,n_want,n_eig,tol,itmax,m_max,use_mat,iseed,both,allsvd
    if (check_lapack) write(10,1200)
    if (check_lapack) write(10,*)
    write(10,*)
    write(10,1300)
    close (10)
    print *
!
!   allocate memory for the eigenvalues and eigenvectors
!
    allocate (eig(n), evec_r(n,n_eig), evec_l(n,n_eig))
!
!   compute a guess for the eigenvector (see guess_evec for more information)
!
  call guess_evec(1,n,n_eig,diagonal,evec_r)
  call dcopy(n*n_eig,evec_r,1,evec_l,1)
!   call driver nonsym
!
  call nonsym_driver(.false.,n,n_want,n_eig,itmax,tol,m_max,0.0d0,mmult,mmult_l,mprec,eig,evec_r,evec_l,both,ok,allsvd)
!
!
  deallocate(eig, evec_r, evec_l, diagonal)
!
  end subroutine test_nonsym
!
  subroutine matexp(n,mat,expmat)
    use real_precision
    implicit none
    integer,                  intent(in)    :: n
    real(dp), dimension(n,n), intent(in)    :: mat 
    real(dp), dimension(n,n), intent(inout) :: expmat
!
    integer  :: i, iter
    real(dp) :: del_norm, fac
    real(dp) :: dnrm2
    external :: dnrm2
!
    real(dp), allocatable :: x(:,:), y(:,:)
!
    integer,  parameter :: maxit = 100
    real(dp), parameter :: tol = 1.0d-20
!
    allocate (x(n,n), y(n,n))
!
    expmat = 0.0_dp
    x      = 0.0_dp
    do i = 1, n
      x(i,i) = 1.0_dp
    end do
!
    do iter = 1, maxit
      del_norm = dnrm2(n*n,x,1)
      expmat = expmat + x
      if (del_norm.lt.tol) exit
      x = matmul(mat,x) / real (iter, kind=dp)
    end do
!
    deallocate (x,y)
    return
  end subroutine matexp
!
  subroutine printMatrix(n,m,A,lda) 
!   
! print formatted matrix
!
    use real_precision
    implicit none
    integer , intent(in)  :: n, m, lda
    real(dp), intent(in)  :: A(lda,lda)
!
    integer :: i, j
!
    do i = 1, n
      do j = 1, m
        write(*,'(F12.5)', advance='no') A(i,j)
        if (j .lt. m) then
          write(*, '(A)', advance='no') ' '
        end if
      end do
      print *
    end do
!
  end subroutine printMatrix
!
  subroutine sort_eigenpairs(wr,wl,vr,vl,n,m,n_want,ldv,ignore,thresh)
!
!   sort m real & imaginary eigenvalues and right & left eigenvectors of length n 
!   in decreasing order according to the real eigenvalues in the range of n_want
! 
    use real_precision
    implicit none
    integer,  intent(in)      :: n, m, ldv, n_want
    real(dp), intent(inout)   :: wr(m), wl(m), vr(ldv,m), vl(ldv,m)
    real(dp), intent(in)      :: thresh
    logical,  intent(in)      :: ignore
!   
!   local variables
!
    real(dp)                  :: w, v(ldv)
    integer                   :: i, j, idx, min_idx(1), fin
    logical                   :: mask(m)
!
    real(dp)                  :: dnrm2
!
    mask = .true.
!
    do i = 1, n_want
! 
!     identify minimal value and mask first position for next iteration
!
print*, wr
      min_idx = minloc(wr, mask=mask) 
      idx     = min_idx(1)
!
!     check complex contribution, if so, move it to with to the last position 
!     of the array and mask it. search again for lowest eigenvalue and 
!     continue with that.
!
      if (ignore .and. abs(wl(idx)) > thresh) then
        fin = m
!
        do j = 1, m
          if (.not. mask(fin)) then
            fin = fin - 1
          else 
            exit
          end if
        end do
!
        mask(fin) = .false.
!
!       do various swaps for double value on last available position fin
!
        w       = wr(fin)
        wr(fin) = wr(idx)
        wr(idx) = w
!
        w       = wl(fin)
        wl(fin) = wl(idx)
        wl(idx) = w
!
        v         = vr(:,fin)
        vr(:,fin) = vr(:,idx)
        vr(:,idx) = v
!    
        v         = vl(:,fin)
        vl(:,fin) = vl(:,idx)
        vl(:,idx) = v
!
!       now search again for lowest and find automatically the corresponding 
!       pair with imaginary contribution
!
        min_idx = minloc(wr, mask=mask) 
        idx     = min_idx(1)
      end if
!
      mask(i) = .false.
!
!     do various swaps to move minimum value et alii on position i
!
      w       = wr(i)
      wr(i)   = wr(idx)
      wr(idx) = w
!
      w       = wl(i)
      wl(i)   = wl(idx)
      wl(idx) = w
!
      v         = vr(:,i)
      vr(:,i)   = vr(:,idx)
      vr(:,idx) = v
!    
      v         = vl(:,i)
      vl(:,i)   = vl(:,idx)
      vl(:,idx) = v
! 
!
    end do
!
  end subroutine sort_eigenpairs
!
  subroutine guess_evec(iwhat,n,m,diagonal,evec)
    use real_precision
    implicit none
    integer,                  intent(in)    :: iwhat, n, m
    real(dp), dimension(n),   intent(in)    :: diagonal
    real(dp), dimension(n,m), intent(inout) :: evec
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
    evec = 0.0_dp
!
    allocate (mask(n))
!
    if (iwhat.eq.1) then
!
!     get the minimum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = minloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = 1.0d0
      enddo
    else if (iwhat.eq.2) then
!
!     get the maximum element of the diagonal
!
      mask = .true.
      do i = 1, m
        ipos = maxloc(diagonal,dim=1,mask=mask)
        mask(ipos) = .false.   
        evec(ipos,i) = 1.0d0
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
        evec(ipos,i) = evec(ipos,i) + 1.0d0
      enddo

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
        evec(ipos,i) = evec(ipos,i) + 1.0d0
      enddo
    end if
    return
  end subroutine guess_evec









