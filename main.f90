program main
  use real_precision
  use utils
  use diaglib
!
! simple minded test program to call the lobpcg routine
!
  implicit none
  integer  :: n, n2, n_eig, itmax, i, j, n_want, m_max, n_seed, ipos
  real(dp) :: shift, tol, lw(1), fac, sqrttwo
  logical  :: ok
  external :: mmult, mprec, smult, apbvec, ambvec, spdvec, smdvec, lrprec_1, lreffprec
!
  integer,  allocatable :: iseed(:)
  real(dp), allocatable :: eig(:), evec(:,:), diagonal(:), a_copy(:,:), s_copy(:,:), w(:)
  real(dp), allocatable :: adiag(:), sigdiag(:) 
  logical,  allocatable :: mask(:)
!
! initialize:
!
  n      = 1000
  n2     = 2*n
  n_want = 6
  n_eig  = min(2*n_want,n_want+5)
  tol    = 1.0e-10_dp
  shift  = 0.0_dp
  itmax  = 1000
  nmult  = 0
  m_max  = 20
!
! allocate memory. note that the matrix a is defined in the stack of module "utils".
!
  allocate (a(n2,n2), eig(n_eig), evec(n2,n_eig), diagonal(n))
  allocate (s(n2,n2))
  call random_seed(size=n_seed)
  allocate (iseed(n_seed))
  iseed = 1
  call random_seed(put=iseed)
  deallocate (iseed)
!
! allocate memory for the a, b, apb, amb, sigma, delta, spd, smd matrices:
!
  allocate (aa(n,n), bb(n,n), apb(n,n), amb(n,n), sigma(n,n), delta(n,n), spd(n,n), smd(n,n))
  allocate (adiag(n), sigdiag(n))
!
! build a positive definite, symmetric apb, amb, sigma matrices:
!
  do i = 1, n
    apb(i,i) = 5.0d0 + dble(i)
    do j = 1, i - 1
      apb(j,i) = 1.0d0/dble(i+j)
      apb(i,j) = apb(j,i)
    end do 
  end do
  do i = 1, n
    amb(i,i) =2.0d0 + dble(i)
    do j = 1, i - 1
      apb(j,i) = 0.2d0/dble(i+j)
      apb(i,j) = apb(j,i)
    end do 
  end do
!
  call random_number(sigma)
  delta = matmul(transpose(sigma),sigma)
  sigma = delta
  do i = 1, n
    sigma(i,i) = sigma(i,i) + 1.0d0
  end do
!
! build antisymmetric delta:
!
  call random_number(delta)
  delta = delta - transpose(delta)
!
! build a and b:
!
  aa = 0.5d0 * (apb + amb)
  bb = 0.5d0 * (apb - amb)
!
! build sigma + delta and sigma - delta
!
  spd = sigma + delta
  smd = sigma - delta 
!
! build the complete matrices:
!
  a(1:n,   1:n)    = aa
  a(n+1:n2,n+1:n2) = aa
  a(1:n,   n+1:n2) = bb
  a(n+1:n2,1:n)    = bb
  s(1:n,   1:n)    = sigma
  s(n+1:n2,n+1:n2) = - sigma
  s(1:n,   n+1:n2) = delta
  s(n+1:n2,1:n)    = - delta
!
! compute the eigenvalues/vectors with lapack for debug
!
  allocate (w(n2), a_copy(n2,n2), s_copy(n2,n2))
  call dsygv(1,'V','L',n2,s_copy,n2,a_copy,n2,w,lw,-1,info)
  lwork = int(lw(1))
  allocate (work(lwork))
  s_copy = s
  a_copy = a
  call dsygv(1,'V','L',n2,s_copy,n2,a_copy,n2,w,work,lwork,info)
  deallocate (work)
 
  open (unit = 10, file = 'lapack.txt', form = 'formatted', access = 'sequential')
  do i = 1, n_want
    write(10,1000) i, 1.0d0/w(n2-i+1)
    if (s_copy(1,n2-i+1) .lt. 0.0d0) s_copy(:,n2-i+1) = - s_copy(:,n2-i+1)
    write(10,'(10f12.6)') s_copy(:,n2-i+1)
    write(10,*)
  end do
 
  close (10)
 
  do i = 1, n
    adiag(i)   = aa(i,i)  
    sigdiag(i) = sigma(i,i)
  enddo
!
  diagonal = adiag - sigdiag
!
  1000 format(t3,' eigenvalue # ',i6,': ',f12.6,/,t3,' eigenvector: ')
!
! test caslr for various guesses.
!
  do i = 2, 6
    call guess_evec(i,n2,n_eig,diagonal,evec)
    call caslr_eff_driver(.true.,n,n2,n_want,n_eig,itmax,tol,m_max,apbvec,ambvec,spdvec,smdvec,lreffprec,eig,evec,ok)
  end do
  open (unit = 10, file = 'newdav.txt', form = 'formatted', access = 'sequential')
  do i = 1, n_want
    write(10,1000) i, eig(i)
    if (evec(1,i) .lt. 0.0d0) evec(:,i) = - evec(:,i)
    write(10,'(10f12.6)') evec(:,i)/2.0d0
    write(10,*)
  end do
!
  close (10)
!
  deallocate (a_copy, s_copy, w, a, s, aa, bb, apb, amb, sigma, delta, spd, smd)
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
        yp(i,icol) = - 1.0d0/(aa(i,i)**2 - fac**2 * sigma(i,i)**2) * (a(i,i) * xp(i,icol) + fac * sigma(i,i) * xm(i,icol))
        ym(i,icol) = - 1.0d0/(aa(i,i)**2 - fac**2 * sigma(i,i)**2) * (a(i,i) * xm(i,icol) + fac * sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lrprec_1
!
  subroutine lreffprec(n,m,fac,xp,xm,yp,ym)
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
        yp(i,icol) = 1.0d0/(fac*fac*aa(i,i)**2 - sigma(i,i)**2) * (fac * a(i,i) * xp(i,icol) + sigma(i,i) * xm(i,icol))
        ym(i,icol) = 1.0d0/(fac*fac*aa(i,i)**2 - sigma(i,i)**2) * (fac * a(i,i) * xm(i,icol) + sigma(i,i) * xp(i,icol))
      end do
    end do
!
    return
  end subroutine lreffprec
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









