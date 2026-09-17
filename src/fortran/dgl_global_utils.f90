module dgl_global_utils
!* Internal module containing constants, the context of a driver call, and the memory
! management, error handling and timing utilities.
!
! DiagLib keeps no global state: everything that belongs to a driver call (memory
! bookkeeping, error status, verbosity) is stored in a dgl_context variable, which is
! local to the driver and passed to the internal routines that need it. This makes
! the drivers re-entrant and thread-safe (as long as the user-supplied routines are).
    use dgl_interface, only: ip => dgl_int, dp => dgl_real, &
                             dgl_success, dgl_err_input, dgl_err_memory, &
                             dgl_err_lapack, dgl_err_ortho, dgl_err_mismatch
    use dgl_lapack
    implicit none
!
    integer, parameter :: i8 = selected_int_kind(18)
!! 64-bit integer kind, always used for memory bookkeeping to avoid overflows
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
!! Useful Constant
    real(dp), parameter :: num_thresh = 1.e-13_dp
!! Numerical threshold for real numbers comparisions with 0
    integer(ip), parameter :: dgl_min_dav_iter = 10_ip
!! Minimum number of iterations between two restarts of the Davidson drivers
!
    type :: dgl_context
!! State of a driver call
        logical :: verbose = .false.
!! Verbosity of the driver call
        integer(ip) :: status = dgl_success
!! Status of the call. Errors do not stop the program: the routine that detects them
!! records the first one here and returns, and the driver cleans up and reports it
!! through its dgl_info argument.
        integer(i8) :: maxmem = 0, maxcor = 0, peakmem = 0
!! Memory bookkeeping (in number of 8-byte words): memory still available,
!! memory available at the beginning, and minimum of the available memory
    end type dgl_context
!
! intrisic functions:
! ===================
!
    intrinsic :: random_number
!
! allocation and deallocations routines
!
!> Overloading for general allocation
    interface mallocate
        module procedure r_alloc1
        module procedure r_alloc2
        module procedure i_alloc1
        module procedure i_alloc2
        module procedure l_alloc1
    end interface mallocate
!
!> Overloading for general deallocation
    interface mfree
        module procedure r_free1
        module procedure r_free2
        module procedure i_free1
        module procedure i_free2
        module procedure l_free1
    end interface mfree
!
contains
!
    subroutine r_alloc1(ctx, len1, v)
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: len1
        real(dp), allocatable, intent(inout) :: v(:)
!
        integer(ip) :: istat
!
        allocate (v(len1), stat=istat)
        call chk_mall(ctx, int(len1, i8), istat)
!
    end subroutine r_alloc1
!
    subroutine r_alloc2(ctx, len1, len2, v)
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: len1, len2
        real(dp), allocatable, intent(inout) :: v(:, :)
!
        integer(ip) :: istat
!
        allocate (v(len1, len2), stat=istat)
        call chk_mall(ctx, int(len1, i8)*int(len2, i8), istat)
!
    end subroutine r_alloc2
!
    subroutine i_alloc1(ctx, len1, v)
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: len1
        integer(ip), allocatable, intent(inout) :: v(:)
!
        integer(ip) :: istat
!
        allocate (v(len1), stat=istat)
        call chk_mall(ctx, int(len1, i8), istat)
!
    end subroutine i_alloc1
!
    subroutine i_alloc2(ctx, len1, len2, v)
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: len1, len2
        integer(ip), allocatable, intent(inout) :: v(:, :)
!
        integer(ip) :: istat
!
        allocate (v(len1, len2), stat=istat)
        call chk_mall(ctx, int(len1, i8)*int(len2, i8), istat)
!
    end subroutine i_alloc2
!
    subroutine l_alloc1(ctx, len1, v)
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: len1
        logical, allocatable, intent(inout) :: v(:)
!
        integer(ip) :: istat
!
        allocate (v(len1), stat=istat)
        call chk_mall(ctx, int(len1, i8), istat)
!
    end subroutine l_alloc1
!
    subroutine r_free1(ctx, v)
        type(dgl_context), intent(inout) :: ctx
        real(dp), allocatable, intent(inout) :: v(:)
!
        integer(i8) :: lfree
        integer(ip) :: istat
!
        if (.not. allocated(v)) return
        lfree = size(v, kind=i8)
        deallocate (v, stat=istat)
        call chk_free(ctx, lfree, istat)
!
    end subroutine r_free1
!
    subroutine r_free2(ctx, v)
        type(dgl_context), intent(inout) :: ctx
        real(dp), allocatable, intent(inout) :: v(:, :)
!
        integer(i8) :: lfree
        integer(ip) :: istat
!
        if (.not. allocated(v)) return
        lfree = size(v, kind=i8)
        deallocate (v, stat=istat)
        call chk_free(ctx, lfree, istat)
!
    end subroutine r_free2
!
    subroutine i_free1(ctx, v)
        type(dgl_context), intent(inout) :: ctx
        integer(ip), allocatable, intent(inout) :: v(:)
!
        integer(i8) :: lfree
        integer(ip) :: istat
!
        if (.not. allocated(v)) return
        lfree = size(v, kind=i8)
        deallocate (v, stat=istat)
        call chk_free(ctx, lfree, istat)
!
    end subroutine i_free1
!
    subroutine i_free2(ctx, v)
        type(dgl_context), intent(inout) :: ctx
        integer(ip), allocatable, intent(inout) :: v(:, :)
!
        integer(i8) :: lfree
        integer(ip) :: istat
!
        if (.not. allocated(v)) return
        lfree = size(v, kind=i8)
        deallocate (v, stat=istat)
        call chk_free(ctx, lfree, istat)
!
    end subroutine i_free2
!
    subroutine l_free1(ctx, v)
        type(dgl_context), intent(inout) :: ctx
        logical, allocatable, intent(inout) :: v(:)
!
        integer(i8) :: lfree
        integer(ip) :: istat
!
        if (.not. allocated(v)) return
        lfree = size(v, kind=i8)
        deallocate (v, stat=istat)
        call chk_free(ctx, lfree, istat)
!
    end subroutine l_free1
!
    subroutine nums_to_bytes(num, b_num, b_unit)
!! Converter from numbers to Bytes. Actual unit depends on the
!! magnitude of the number and it is returned. Assumes 8 Byte numbers
        implicit none
        integer(i8), intent(in) :: num
        real(dp), intent(inout) :: b_num
        character(len=*), intent(inout) :: b_unit
!
        integer(i8) :: num_l
        real(dp) :: converter
!
        num_l = 8_i8*num
        select case (num_l)
        case (:1000000_i8)
            b_unit = "KB"
            converter = 1.e3_dp
!
        case (1000001_i8:1000000000_i8)
            converter = 1.e6_dp
            b_unit = "MB"
!
        case default
            converter = 1.e9_dp
            b_unit = "GB"
        end select
!
        b_num = real(num_l, kind=dp)/converter
!
    end subroutine nums_to_bytes
!
    subroutine bytes_to_nums(verbose, bytes, bytes_unit, nums)
!! Converter from Bytes to numbers. Assumes 8 Byte numbers
        implicit none
        logical, intent(in) :: verbose
        integer(ip), intent(in) :: bytes
        character(len=*), intent(in) :: bytes_unit
        integer(i8), intent(inout) :: nums
!
        integer(i8) :: converter
!
        select case (bytes_unit)
        case ("KB")
            converter = 1000_i8
        case ("MB")
            converter = 1000000_i8
        case ("GB")
            converter = 1000000000_i8
        case default
            if (verbose) write (*, "(t3,a,/)") "Unknown or Unspecified memory_unit, defaulting to MBs"
            converter = 1000000_i8
        end select
!
        nums = int(bytes, i8)*converter/8_i8
!
    end subroutine bytes_to_nums
!
    subroutine chk_mall(ctx, lall, istat)
!! Check for proper allocations. Also keeps track of the memory used.
!! Checks for out of memory condition.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(i8), intent(in) :: lall
        integer(ip), intent(in) :: istat
!
        real(dp) :: b_lall, b_maxmem
        character(len=2) :: lall_unit, maxmem_unit
!
9000    format(t3, 'allocation error, stat= ', i5)
9010    format(t3, 'allocation error,', /, &
               t3, 'not enough memory. ', f10.3, ' ', a2, ' required', /, &
               t3, '                   ', f10.3, ' ', a2, ' available.')
!
        if (istat .ne. 0) then
            if (.not. dgl_failed(ctx)) write (*, 9000) istat
            call dgl_error(ctx, "allocation failed", dgl_err_memory)
            return
        end if
!
! the array has been allocated: account for it even when the memory limit is
! exceeded, so that the bookkeeping stays consistent when it is freed.
!
        if (lall .gt. ctx%maxmem) then
            if (.not. dgl_failed(ctx)) then
                call nums_to_bytes(lall, b_lall, lall_unit)
                call nums_to_bytes(ctx%maxmem, b_maxmem, maxmem_unit)
                write (*, 9010) b_lall, lall_unit, b_maxmem, maxmem_unit
            end if
            call dgl_error(ctx, "memory limit exceeded", dgl_err_memory)
        end if
        ctx%maxmem = ctx%maxmem - lall
        if (ctx%peakmem .gt. ctx%maxmem) ctx%peakmem = ctx%maxmem
!
    end subroutine chk_mall
!
    subroutine chk_free(ctx, lfree, istat)
!! Check for proper deallocation.
!! Also keeps track of memory released.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(i8), intent(in) :: lfree
        integer(ip), intent(in) :: istat
!
9000    format(t3, 'deallocation error, stat= ', i5)
!
        if (istat .ne. 0) then
            if (.not. dgl_failed(ctx)) write (*, 9000) istat
            call dgl_error(ctx, "deallocation failed", dgl_err_memory)
        else
            ctx%maxmem = ctx%maxmem + lfree
        end if
    end subroutine chk_free
!
    subroutine dgl_check_memleak(ctx)
!! Check for the presence of internal memory leaks
        implicit none
        type(dgl_context), intent(in) :: ctx
        real(dp) :: p_mem
        character(len=2) :: p_mem_unit
!
! check if DiagLib is leaking data
!
        if (ctx%maxmem .ne. ctx%maxcor) then
            write (*, "(t3,a,/,t3,i0,a,/)") "DiagLib: Memory leak detected. DiagLib should free ", &
                ctx%maxcor - ctx%maxmem, " more numbers from memory"
        end if
        call nums_to_bytes(ctx%maxcor - ctx%peakmem, p_mem, p_mem_unit)
        if (ctx%verbose) write (*, "(t3,a,f10.3,a3,/)") "DiagLib peak memory used:", p_mem, p_mem_unit
!
    end subroutine dgl_check_memleak
!
! =====================
! More global routines
! =====================
!
    subroutine dgl_init(ctx, lenght, n_arrs, mem, mem_unit, verbose_in)
!! Initializer for all drivers: sets the verbosity and the memory available to the call
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: lenght, n_arrs
        integer(ip), intent(in) :: mem
        character(len=2), intent(in) :: mem_unit
        logical, intent(in) :: verbose_in
!
        integer(i8) :: numbers, numbers_preview
        real(dp) :: memory_preview
        character(len=2) :: memory_preview_unit
!
        ctx%verbose = verbose_in
!
! set maximum memory used by DiagLib to the input value
!
        numbers_preview = int(lenght, i8)*int(n_arrs, i8)
        call nums_to_bytes(numbers_preview, memory_preview, memory_preview_unit)
        if (ctx%verbose) write (*, "(t3,a,f10.3,a3,/)") "DiagLib estimated memory usage is", &
            memory_preview, memory_preview_unit
!
        call bytes_to_nums(ctx%verbose, mem, mem_unit, numbers)
        if (numbers .lt. numbers_preview) call dgl_warning("provided memory is probably not sufficient")
        ctx%maxcor = numbers
        ctx%maxmem = ctx%maxcor
        ctx%peakmem = ctx%maxcor
!
    end subroutine dgl_init
!
    subroutine dgl_error(ctx, string, code)
!! Record a DiagLib error. Only the first error of a driver call is printed and kept;
!! the caller is responsible for returning, and the driver for reporting it.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        character(len=*), intent(in) :: string
        integer(ip), intent(in) :: code

        if (ctx%status .ne. dgl_success) return
        ctx%status = code
        write (*, "(t3,a)") "-- DiagLib Error: "//string
    end subroutine
!
    subroutine dgl_check_input(ctx, n, n_targ, n_max, max_iter, tol, memory)
!! Checks on the input common to all the drivers. Errors are recorded with dgl_error.
        implicit none
        type(dgl_context), intent(inout) :: ctx
        integer(ip), intent(in) :: n
!! Size of the problem
        integer(ip), intent(in) :: n_targ, n_max, max_iter, memory
        real(dp), intent(in) :: tol
!
        if (n .lt. 1) call dgl_error(ctx, "The size of the problem must be positive", dgl_err_input)
        if (n_targ .lt. 1) call dgl_error(ctx, "The number of requested eigenpairs (n_targ) must be positive", &
                                          dgl_err_input)
        if (n_targ .gt. n_max) call dgl_error(ctx, &
                                              "Number of eigenvalues requested is larger that size of arrays passed", &
                                              dgl_err_input)
        if (max_iter .lt. 1) call dgl_error(ctx, "The maximum number of iterations must be positive", dgl_err_input)
!
! written this way to catch a NaN as well
!
        if (.not. tol .gt. zero) call dgl_error(ctx, "The convergence threshold must be positive", dgl_err_input)
        if (memory .lt. 1) call dgl_error(ctx, "The memory available to DiagLib must be positive", dgl_err_input)
!
! the number of elements of the n x n_max arrays is passed to blas routines as an integer
!
        if (int(n, i8)*int(n_max, i8) .gt. int(huge(n), i8)) call dgl_error(ctx, &
            "Problem too large for the integers of this DiagLib build: use a build with 64-bit integers", &
            dgl_err_input)
!
    end subroutine dgl_check_input
!
    logical function dgl_failed(ctx)
!! True if an error has been recorded in the driver call
        implicit none
        type(dgl_context), intent(in) :: ctx
        dgl_failed = ctx%status .ne. dgl_success
    end function dgl_failed
!
    subroutine dgl_return_info(ctx, info)
!! Report the error status of a driver call: through info, if present,
!! otherwise by stopping the program on error.
        implicit none
        type(dgl_context), intent(in) :: ctx
        integer(ip), optional, intent(out) :: info

        if (present(info)) then
            info = ctx%status
        else if (ctx%status .ne. dgl_success) then
            write (*, "(t3,a)") "** DiagLib issued stop signal **"
            stop 1
        end if
    end subroutine dgl_return_info
!
    subroutine dgl_warning(string)
!! Print a DiagLib warning
        implicit none
        character(len=*), intent(in) :: string

        write (*, "(t3,a)") "-- DiagLib Warning: "//string
    end subroutine
!
    integer(ip) function get_mem_lapack(len_red)
!! Size of the Lapack work array needed by the drivers to diagonalize reduced
!! matrices of dimension up to len_red, with dsyev (symmetric) or dgeev
!! (non-symmetric, with eigenvectors).
!! The minimum requirements are 3*len_red-1 and 4*len_red, respectively;
!! the block sizes returned by ilaenv allow for the blocked algorithms.
        implicit none
        integer(ip), intent(in) :: len_red
!! Maximum dimension of the reduced matrices
        integer(ip) :: nb
!
        nb = max(1_ip, ilaenv(1_ip, 'DSYTRD', 'U', len_red, -1_ip, -1_ip, -1_ip), &
                 ilaenv(1_ip, 'DGEHRD', ' ', len_red, 1_ip, len_red, -1_ip), &
                 ilaenv(1_ip, 'DORGHR', ' ', len_red, 1_ip, len_red, -1_ip))
!
        get_mem_lapack = max(1_ip, (nb + 2_ip)*len_red, 4_ip*len_red + 2_ip*nb*len_red)
    end function get_mem_lapack
!
    subroutine get_time(t)
!! Stores cpu/wall time
!$      use omp_lib, only: omp_get_wtime
        real(dp), dimension(2), intent(inout) :: t
!
! get cpu and (if openmp is available) wall time.
!
        t = zero
        call cpu_time(t(1))
!$      t(2) = omp_get_wtime()
!
        return
    end subroutine get_time

end module dgl_global_utils
