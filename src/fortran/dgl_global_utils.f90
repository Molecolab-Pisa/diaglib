module dgl_global_utils
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
! variables to keep track of memory
!
  integer, private :: maxmem, maxcor
!
! external functions:
! ===================
!
  real(dp) :: dnrm2
  external dgeqrf, dtrsm, dpotrf, dtrcon, dgemm, dnrm2
!
! intrisic functions:
! ===================
!
  intrinsic :: random_number
!
! allocation and deallocations routines
!
  interface mallocate
    module procedure r_alloc1
    module procedure r_alloc2
    module procedure i_alloc1
    module procedure i_alloc2
    module procedure c_alloc1
    module procedure c_alloc2
    module procedure ch_alloc1
    module procedure l_alloc1
  end interface mallocate
!
  interface mfree
    module procedure r_free1
    module procedure r_free2
    module procedure i_free1
    module procedure i_free2
    module procedure c_free1
    module procedure c_free2
    module procedure ch_free1
    module procedure l_free1
  end interface mfree
!
  contains
!
  subroutine r_alloc1(len1,v)
    implicit none
    integer,               intent(in)    :: len1
    real(dp), allocatable, intent(inout) :: v(:)
!
    integer :: istat
!    
    allocate (v(len1), stat=istat)
    call chk_mall(len1,istat)
!
  end subroutine r_alloc1
!
  subroutine r_alloc2(len1,len2,v)
    implicit none
    integer,               intent(in)    :: len1, len2
    real(dp), allocatable, intent(inout) :: v(:,:)
!
    integer :: istat
!    
    allocate (v(len1,len2), stat=istat)
    call chk_mall(len1*len2,istat)
!
  end subroutine r_alloc2
!
  subroutine i_alloc1(len1,v)
    implicit none
    integer,              intent(in)    :: len1
    integer, allocatable, intent(inout) :: v(:)
!
    integer :: istat
!    
    allocate (v(len1), stat=istat)
    call chk_mall(len1,istat)
!
  end subroutine i_alloc1
!
  subroutine i_alloc2(len1,len2,v)
    implicit none
    integer,              intent(in)    :: len1, len2
    integer, allocatable, intent(inout) :: v(:,:)
!
    integer :: istat
!    
    allocate (v(len1,len2), stat=istat)
    call chk_mall(len1*len2,istat)
!
  end subroutine i_alloc2
!
  subroutine c_alloc1(len1,v)
    implicit none
    integer,                  intent(in)    :: len1
    complex(dp), allocatable, intent(inout) :: v(:)
!
    integer :: istat
!    
    allocate (v(len1), stat=istat)
    call chk_mall(len1,istat)
!
  end subroutine c_alloc1
!
  subroutine c_alloc2(len1,len2,v)
    implicit none
    integer,                  intent(in)    :: len1, len2
    complex(dp), allocatable, intent(inout) :: v(:,:)
!
    integer :: istat
!    
    allocate (v(len1,len2), stat=istat)
    call chk_mall(len1*len2,istat)
!
  end subroutine c_alloc2
!
  subroutine ch_alloc1(len1,v)
    implicit none
    integer,                       intent(in)    :: len1
    character(len=*), allocatable, intent(inout) :: v(:)
!
    integer :: istat
!    
    allocate (v(len1), stat=istat)
    call chk_mall(len1,istat)
!
  end subroutine ch_alloc1
!
  subroutine l_alloc1(len1,v)
    implicit none
    integer,              intent(in)    :: len1
    logical, allocatable, intent(inout) :: v(:)
!
    integer :: istat
!    
    allocate (v(len1), stat=istat)
    call chk_mall(len1,istat)
!
  end subroutine l_alloc1
!
  subroutine r_free1(v)
    real(dp), allocatable, intent(inout) :: v(:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine r_free1
!
  subroutine r_free2(v)
    real(dp), allocatable, intent(inout) :: v(:,:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine r_free2
!
  subroutine i_free1(v)
    integer, allocatable, intent(inout) :: v(:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine i_free1
!
  subroutine i_free2(v)
    integer, allocatable, intent(inout) :: v(:,:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine i_free2
!
  subroutine c_free1(v)
    complex(dp), allocatable, intent(inout) :: v(:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine c_free1
!
  subroutine c_free2(v)
    complex(dp), allocatable, intent(inout) :: v(:,:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine c_free2
!
  subroutine ch_free1(v)
    character(len=*), allocatable, intent(inout) :: v(:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine ch_free1
!
  subroutine l_free1(v)
    logical, allocatable, intent(inout) :: v(:)
!
    integer :: lfree, istat
!
    if (.not. allocated(v)) return
    lfree = size(v)
    deallocate (v, stat=istat)
    call chk_free(lfree,istat)
!
  end subroutine l_free1
!
!
  subroutine to_xbytes(num,b_num,b_unit)
    implicit none
    integer,      intent(in) :: num
    real(dp),         intent(inout) :: b_num
    character(len=*), intent(inout) :: b_unit
!
    integer :: num_l
!
    num_l = 8*num
    select case(num_l)
    case(:1000000)
      b_num = real(num,kind=dp) / 1.e3_dp
      b_unit = "KB"
    case (1000001:1000000000)
      b_num = real(num,kind=dp) / 1.e6_dp
      b_unit = "MB"
    case default
      b_num = real(num,kind=dp) / 1.e9_dp
      b_unit = "GB"
    end select
!
  end subroutine to_xbytes
!
  subroutine chk_mall(lall,istat)
    implicit none
    integer,           intent(in) :: lall,  istat
!
    real(dp) :: b_lall, b_maxmem
    character(len=2) :: lall_unit, maxmem_unit
!
9000 format(t3,'allocation error, stat= ',i5)
9010 format(t3,'allocation error,',/,    &
            t3,'not enough memory. ',f10.3,' ',a2,' required',/, &
            t3,'                   ',f10.3,' ',a2,' available.')
!
    if (istat.ne.0) then
      write(*,9000) istat
      stop
    else if (lall.gt.maxmem) then
      call to_xbytes(lall,b_lall,lall_unit)
      call to_xbytes(maxmem,b_maxmem,maxmem_unit)
      write(*,9010) b_lall, lall_unit, b_maxmem, maxmem_unit
      stop
    else
      maxmem = maxmem - lall
    end if
!    
  end subroutine chk_mall
!
  subroutine chk_free(lfree,istat)
    implicit none
    integer,           intent(in) :: lfree, istat
!
 9000 format(t3,'deallocation error, stat= ',i5)
!
    if (istat.ne.0) then
      write(*,9000) istat
      stop
    else
      maxmem = maxmem + lfree
    end if
  end subroutine chk_free
!
  subroutine dgl_init_(mem)
    implicit none
    integer, intent(in) :: mem
!
! set maximum memory used by DiagLib to the input value
!
    maxcor = mem
    maxmem = maxcor
!
  end subroutine dgl_init_
!
  subroutine dgl_check_memleak()
    implicit none
!
! check if DiagLib is globally leaking data
!
    if (maxmem.ne.maxcor) then
      write(*,"(t3,a,/,t3,i0,a,/)") "Memory leak detected. DiagLib should free ", maxcor-maxmem, &
                    " more words of space"
    endif
!
  end subroutine dgl_check_memleak
!
!  =====================
!  More global routines
!  =====================
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

end module dgl_global_utils
