module dgl_global_utils
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
!! Global variable holding kind for double precision
  real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
!! Useful Constant
  integer :: lwork, info
!! Lapack utility  
  real(dp), parameter   :: num_thresh = 1.e-13_dp
!! Numerical threshold for real numbers comparisions with 0
  real(dp), allocatable  :: work(:), tau(:)
!! Lapack utility
  real(dp)               :: t1(2), t2(2), t_diag(2), t_ortho(2), &
                            t_mv(2), t_tot1(2), t_tot2(2), t_tot(2)
!! Timings
  integer, protected :: maxmem, maxcor, peakmem
!! Variables to keep track of memory
  logical, protected :: verbose
!! Global verbosity mode
!
! external functions:
! ===================
!
  real(dp) :: dnrm2
!! Lapack  
  external dgeqrf, dtrsm, dpotrf, dtrcon, dgemm, dnrm2
!! Lapack  
!
! intrisic functions:
! ===================
!
  intrinsic :: random_number
!! Intrinsic function  
!
! allocation and deallocations routines
!
!> Overloading for general allocation
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
!> Overloading for general deallocation
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
  subroutine nums_to_bytes(num,b_num,b_unit)
  !! Converter from numbers to Bytes. Actual unit depends on the
  !! magnitude of the number and it is returned. Assumes 8 Byte numbers
    implicit none
    integer,          intent(in) :: num
    real(dp),         intent(inout) :: b_num
    character(len=*), intent(inout) :: b_unit
!
    integer  :: num_l
    real(dp) ::converter
!
    num_l = 8*num
    select case(num_l)
    case(:int(1.e6_dp))
      b_unit = "KB"
      converter = 1.e3_dp
!
    case (int(1.e6_dp)+1:int(1.e9_dp))
      converter = 1.e6_dp
      b_unit = "MB"
!
    case default
      converter = 1.e9_dp
      b_unit = "GB"
    end select
!    
      b_num = real(num_l,kind=dp) / converter
!
  end subroutine nums_to_bytes
!
  subroutine bytes_to_nums(bytes,bytes_unit,nums)
  !! Converter from Bytes to numbers. Assumes 8 Byte numbers
    implicit none
    integer,          intent(in) :: bytes
    character(len=*), intent(in) :: bytes_unit
    integer,          intent(inout) :: nums
!
    integer :: converter
!
    select case(bytes_unit)
    case("KB")
      converter = int(1.e3_dp)
    case ("MB")
      converter = int(1.e6_dp)
    case ("GB")
      converter = int(1.e9_dp)
    case default
      if(verbose) write(*,"(t3,a,/)") "Unknown or Unspecified memory_unit, defaulting to MBs"
      converter = int(1.e6_dp)
    end select
!
    nums = int(bytes * converter / 8)
!
  end subroutine bytes_to_nums
!
  subroutine chk_mall(lall,istat)
!! Check for proper allocations. Also keeps track of the memory used.
!! Checks for out of memory condition.  
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
      call nums_to_bytes(lall,b_lall,lall_unit)
      call nums_to_bytes(maxmem,b_maxmem,maxmem_unit)
      write(*,9010) b_lall, lall_unit, b_maxmem, maxmem_unit
      stop
    else
      maxmem = maxmem - lall
      if (peakmem .gt. maxmem) peakmem = maxmem 
    end if
!    
  end subroutine chk_mall
!
  subroutine chk_free(lfree,istat)
!! Check for proper deallocation.
!! Also keeps track of memory released.
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
  subroutine dgl_check_memleak()
!! Check for the presence of internal memory leaks
    implicit none
    real(dp) :: p_mem
    character(len=2) :: p_mem_unit
!
! check if DiagLib is globally leaking data
!
    if (maxmem.ne.maxcor) then
      write(*,"(t3,a,/,t3,i0,a,/)") "DiagLib: Memory leak detected. DiagLib should free ", maxcor-maxmem, &
                    " more numbers from memory"
    endif
    call nums_to_bytes(maxcor - peakmem,p_mem,p_mem_unit)
    if (verbose) write(*,"(t3,a,f10.3,a3,/)") "DiagLib peak memory used:", p_mem, p_mem_unit
!
  end subroutine dgl_check_memleak
!
!  =====================
!  More global routines
!  =====================
!
  subroutine dgl_init(lenght,n_arrs,mem,mem_unit,verbose_in)
!! Global initializer for all drivers  
    implicit none
    integer,          intent(in) :: lenght, n_arrs
    integer,          intent(in) :: mem
    character(len=2), intent(in) :: mem_unit
    logical,          intent(in) :: verbose_in
!
    integer :: numbers, numbers_preview
    real(dp) :: memory_preview
    character(len=2) :: memory_preview_unit 
!
! Set global verbosity level to user input
!
    verbose = verbose_in
!
! set maximum memory used by DiagLib to the input value
!
    numbers_preview = lenght*n_arrs
    call nums_to_bytes(numbers_preview,memory_preview,memory_preview_unit)
    if (verbose) write(*,"(t3,a,f10.3,a3,/)") "DiagLib esitmated memory usage is", memory_preview, memory_preview_unit
!
    call bytes_to_nums(mem,mem_unit,numbers)
    if (numbers .lt. numbers_preview) write(*,"(t3,a)") "-- DiagLib Warning: provided memory is probably not sufficient"
    maxcor = numbers
    maxmem = maxcor
    peakmem = maxcor
!
! Reset timings
!
    t_tot       = zero
    t_diag      = zero
    t_ortho     = zero
    t_mv        = zero
!
  end subroutine dgl_init
!
  subroutine dgl_error(string)
    implicit none
    character(len=*), intent(in) :: string

    write(*,"(t3,a)") string
    stop "DiagLib isseued stop signal"
  end subroutine
!
  integer function get_mem_lapack(n,n_max)
!! Get the highest optimal memory amount required by lapacks
!! for the driver execution.
    integer, intent(in)    :: n
!! Number of rows of matrices that will be processed
    integer, intent(in)    :: n_max
!! Number of columns of matrices that will be processed
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
!! Stores cpu/wall time  
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
