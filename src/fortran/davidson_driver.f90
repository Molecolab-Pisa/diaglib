subroutine davidson_driver(n,n_targ,n_max,matvec,precnd,eig,evec,ok, &
              dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
              dgl_shift, dgl_memory,dgl_memory_unit, metvec)
!! ### Driver for Davidson-Liu symmetric diagonalization
!! Can solve both standard and generalized eigenvalue problems.
!! In the latter case you need to pass the optional argument [[metvec]] as a pointer to your routine.
!! Moreover, you need to use the [[dgl_drivers_interfaces]] module included in this library.
!!
!! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
  use dgl_minor_utils
  use dgl_external_interfaces
  implicit none
    integer,                      intent(in)    :: n
!! Size of the matrix to be diagonalized
    integer,                      intent(in)    :: n_targ
!! Number of required eigenpairs.
    integer,                      intent(in)    :: n_max
!! Maximum size of the search space. Should be >= n_targ.
    real(dp), dimension(n_max),   intent(inout) :: eig
!! Computed eigenvalues    
    real(dp), dimension(n,n_max), intent(inout) :: evec
!! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
    logical,                      intent(inout) :: ok
!! True if davidson converged
    procedure(matvec_) :: matvec
!! External subroutine that performs the matrix-vector multiplication
    procedure(precnd_) :: precnd
!! External subroutine that applies a preconditioner
    logical,  optional,            intent(in)    :: dgl_verbose
!! Verbose mode. Default = .false.
    integer,  optional,            intent(in)    :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
    integer,  optional,            intent(in)    :: dgl_dav_iter
!! Maximum number of iterations before Davidson restart. Default = \(25\)
    integer,  optional,            intent(in)    :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
    character(len=2),  optional,   intent(in)    :: dgl_memory_unit
!! Unit of memory. Default = MBs
    real(dp), optional,            intent(in)    :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
    real(dp), optional,            intent(in)    :: dgl_shift
!! Diagonal level shifting parameter. Default = \(0.\)
    procedure(metvec_), pointer, optional :: metvec
!! Pointer to External subroutine that applies the metric-vector multiplication
!
!   local variables:
!   ================
    logical  :: verbose_l
    integer  :: max_iter, dav_iter, memory
    real(dp) :: tol, shift
    character(len=2) :: memory_unit
!
!   expansion space varibles: 
!     total dimension, current dimension
!
    integer               :: lda, ld_current
!
!   number of large arrays that will be allocated
!
   integer               :: n_arrs
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
!   tolerances on residuals norms, used for convergence
!
    real (dp)             :: tol_rms, tol_max
!
!   variables for the restart
!
    integer               :: n_rst
    logical               :: restart
!
!   type of problem (standard or generalized)
!    
    logical               :: generalized
!
!   iterators and utilities
!
    integer               :: it, i_eig
    real(dp)              :: sqrtn, xx(1)
!
!   array to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: space(:,:), aspace(:,:), residuals(:,:), r_norm(:,:)
    real(dp), allocatable :: bspace(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), e_red(:)
    real(dp), allocatable :: s_red(:,:), s_copy(:,:), b_evec(:,:)
!
!   ================
!   START EXECUTION
!   ================
!
!  Stupidity check
!
    if(n_targ.gt.n_max) call dgl_error(&
    "Number of eigenvalues request is larger that size of arrays passed")
!
!   check what problem we are dealing with
!
    generalized = present(metvec)
    if (generalized) then
      if (.not.associated(metvec)) then
        stop "DiagLib: Non associated pointer to metric-vector product routine"
      endif
    endif
!
! Parse optional arguments
!
    verbose_l = .false. ; if(present(dgl_verbose)) verbose_l = dgl_verbose
    max_iter = 100      ; if(present(dgl_max_iter)) max_iter = dgl_max_iter
    dav_iter = 25       ; if(present(dgl_dav_iter)) dav_iter = dgl_dav_iter
    tol = 1.e-7_dp      ; if(present(dgl_tol)) tol = dgl_tol
    shift = 0.e0_dp     ; if(present(dgl_shift)) shift = dgl_shift
    memory = 80         ; if(present(dgl_memory)) memory = dgl_memory
    memory_unit = "MB"  ; if(present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
!   compute the actual size of the expansion space
!
    lda = dav_iter*n_max
!
    if(generalized) then
      n_arrs = lda*3 + n_max*2
    else
      n_arrs = lda*2 + n_max
    endif
    call dgl_init(n,n_arrs,memory,memory_unit,verbose_l)
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    call mallocate(lwork,work)
    call mallocate(n_max,tau)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residuals:
!
    call mallocate(n,lda,space)
    call mallocate(n,lda,aspace)
    call mallocate(n,n_max,residuals)

    if (generalized) then
      call mallocate(n,lda,bspace)
      call mallocate(n,n_max,b_evec)
    endif
!
!   allocate memory for convergence check
!
    call mallocate(n_max,done)
    call mallocate(2,n_max,r_norm)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    call mallocate(lda,lda,a_red)
    call mallocate(lda,lda,a_copy)
    call mallocate(lda,e_red)
    if (generalized) then
      call mallocate(lda,lda,s_red)
      call mallocate(lda,lda,s_copy)
    endif
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag   = zero
    t_ortho  = zero
    t_mv     = zero
    t_tot    = zero
    space   = zero
    aspace  = zero
    a_red   = zero
    if (generalized) then
      bspace  = zero
      s_red   = zero
    endif
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess(n,n_max,evec)
!
!   move the guess into the expansion space.
!
    call dcopy(n*n_max,evec,1,space,1)
!
!   apply the b matrix and b-orthogonalize the guess:
!
    if (generalized) then
      call metvec(n,n_max,space,bspace)
      call b_ortho(n,n_max,space,bspace)
    endif
!
!   initialize the number of active vectors and the associated indices.
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    ld_current   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Davidson-Liu iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ld_current = ld_current + n_act
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call matvec(n,n_act,space(1,i_beg+n_rst),aspace(1,i_beg+n_rst))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ld_current,n_act,n,one,space,n,aspace(1,i_beg+n_rst),n,zero,a_red(1,i_beg+n_rst),lda)
!
!     explicitly putting the first block of 
!     converged eigenvalues in the reduced matrix
!
      if (restart) then
        do i_eig = 1, n_rst
          a_red(i_eig,i_eig) = e_red(i_eig)
        end do
        restart = .false.
        n_rst   = 0
      end if
      a_copy = a_red
!
!     diagonalize the reduced matrix
!
      call get_time(t1)
      call dsyev('v','u',ld_current,a_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      eig = e_red(1:n_max)
!
      call dgemm('n','n',n,n_max,ld_current,one,space,n,a_copy,lda,zero,evec,n)
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ld_current,one,aspace,n,a_copy,lda,zero,residuals,n)
      if (generalized) call dgemm('n','n',n,n_max,ld_current,one,bspace,n,a_copy,lda,zero,b_evec,n)
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        if (generalized) then
          call daxpy(n,-eig(i_eig),b_evec(:,i_eig),1,residuals(:,i_eig),1)
        else
          call daxpy(n,-eig(i_eig),evec(:,i_eig),1,residuals(:,i_eig),1)
        endif
        r_norm(1,i_eig) = dnrm2(n,residuals(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(residuals(:,i_eig)))
      end do
!
!     check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*) 
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        exit
      end if
!
!     check whether an update is required. 
!     if not, perform a davidson restart.
!
      if (ld_current .lt. lda) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        i_beg = i_beg + n_act
        n_act = n_max
        n_frozen = 0
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_act = n_act - 1
            n_frozen = n_frozen + 1
          else
            exit
          end if
        end do
        ind   = n_max - n_act + 1
        call precnd(n,n_act,-eig(ind),residuals(1,ind),space(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        if (generalized) then
          call b_ortho_vs_x(n,ld_current,n_act,space,bspace,space(1,i_beg))
          call metvec(n,n_act,space(1,i_beg),bspace(1,i_beg))
          call b_ortho(n,n_act,space(1,i_beg),bspace(1,i_beg))
        else
            call ortho_vs_x(n,ld_current,n_act,space,space(1,i_beg),xx,xx)
        endif
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a,/)') 'Restarting davidson.'
        n_act = n_max
        space = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        call dcopy(n_max*n,evec,1,space,1)
        aspace = zero
        a_red  = zero
        if (generalized) then
          call dcopy(n_max*n,b_evec,1,bspace,1)
          call b_ortho(n,n_max,space,bspace)
          bspace = zero
          s_red  = zero
        endif
!
!       initialize indexes back to their starting values 
!
        ld_current   = 0
        i_beg = 1 
        n_rst = 0
!
!       counting how many matvec we can skip at the next
!       iteration
!
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_rst = n_rst + 1
          else
            exit
          end if
        end do
        restart = .true.
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!
    end do
!
!   deallocate memory
!
    call mfree(work)
    call mfree(tau)
    call mfree(space)
    call mfree(aspace)
    call mfree(bspace)
    call mfree(residuals)
    call mfree(done)
    call mfree(r_norm)
    call mfree(a_red)
    call mfree(a_copy)
    call mfree(e_red)
    if (generalized) then
      call mfree(bspace)
      call mfree(b_evec)
      call mfree(s_red)
      call mfree(s_copy)
    endif
!
    call dgl_check_memleak()
!
!   if required, print timings
!
    call get_time(t2)
    t_tot = t2 - t_tot
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!
1000 format(t3,'timings for Davidson-Liu (cpu/wall): ',/, &
            t3,'  matrix-vector multiplications: ',2f12.4,/, &
            t3,'  diagonalization:               ',2f12.4,/, &
            t3,'  orthogonalization:             ',2f12.4,/, &
            t3,'                                 ',24('='),/,  &
            t3,'  total:                         ',2f12.4)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
!
end subroutine davidson_driver