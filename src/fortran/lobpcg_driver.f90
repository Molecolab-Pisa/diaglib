module mod_lobpcg_driver
  use dgl_minor_utils
  use dgl_external_interfaces

implicit none

contains

subroutine lobpcg_driver(n,n_targ,n_max,matvec,precnd,eig,evec,ok, &
              dgl_verbose, dgl_max_iter, dgl_tol, &
              dgl_shift, dgl_memory,dgl_memory_unit, metvec)
!! ### Driver for LOBPCG symmetric diagonalization
!! Can solve both standard and generalized eigenvalue problems.
!! In the latter case you need to pass the optional argument [[metvec]] as a pointer to your routine.
!! @note
!! eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
!! @endnote
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
    integer,  optional,            intent(in)    :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
    character(len=2),  optional,   intent(in)    :: dgl_memory_unit
!! Unit of memory. Default = MB
    real(dp), optional,            intent(in)    :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
    real(dp), optional,            intent(in)    :: dgl_shift
!! Diagonal level shifting parameter. Default = \(0.\)
    procedure(metvec_), pointer, optional :: metvec
!! Pointer to External subroutine that applies the metric-vector multiplication
!
!   local variables:
!   ================
    logical  :: bool
    integer  :: max_iter, memory
    real(dp) :: tol, shift
    character(len=2) :: memory_unit
!
!   expansion space varibles: 
!       total dimension, current dimension
!
    integer               :: lda, ld_current
!
!   number of large arrays that will be allocated
!
   integer               :: n_arrs
!
!   tolerances on residuals norms, used for convergence
!
    real (dp)             :: tol_rms, tol_max
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act
!
!   indexes to access specific parts of the expansion space
!
    integer               :: ind_x, ind_w, ind_p
!
!   varible to determine the type of problem
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
    real(dp), allocatable :: a_red(:,:), e_red(:)
    real(dp), allocatable :: u_x(:,:), u_p(:,:), x_new(:,:), ax_new(:,:)
    real(dp), allocatable :: bx_new(:,:)
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
    bool = .false.     ; if(present(dgl_verbose)) bool = dgl_verbose
    max_iter = 50      ; if(present(dgl_max_iter)) max_iter = dgl_max_iter
    tol = 1.e-7_dp     ; if(present(dgl_tol)) tol = dgl_tol
    shift = 0.e0_dp    ; if(present(dgl_shift)) shift = dgl_shift
    memory= 80         ; if(present(dgl_memory)) memory = dgl_memory !80MBs
    memory_unit = "MB" ; if(present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
!   set size of the expansion space
!
    lda = 3*n_max
!
    if(generalized) then
      n_arrs = lda*3 + n_max*4
    else
      n_arrs = lda*2 + n_max*3
    endif
    call dgl_init(n,n_arrs,memory,memory_unit,bool)
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    call mallocate(lwork,work)
    call mallocate(2*n_max,tau)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residuals:
!
    call mallocate(n,lda,space)
    call mallocate(n,lda,aspace)
    call mallocate(n,n_max,residuals)
    if (generalized) call mallocate(n,lda,bspace)

!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    call mallocate(lda,lda,a_red)
    call mallocate(lda,e_red)
!
!   allocate memory for temporary copies of x, ax, and bx:
!
    call mallocate(n,n_max,x_new)
    call mallocate(n,n_max,ax_new)
    if (generalized) call mallocate(n,n_max,bx_new)
!
!   allocate memory for convergence check
!
    call mallocate(n_max,done)
    call mallocate(2,n_max,r_norm)
!
!   clean out:
!
    space    = zero
    aspace   = zero
    a_red    = zero
    if (generalized) bspace = zero
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess(n,n_max,evec)
!
!   if required, compute b*evec and b-orthogonalize the guess
!
    if (generalized) then 
      call metvec(n,n_max,evec,bx_new)
      call b_ortho(n,n_max,evec,bx_new)
    end if
!
!   compute the first eigenpairs by diagonalizing the reduced matrix:
!
    call dcopy(n*n_max,evec,1,space,1)
    if (generalized) call dcopy(n*n_max,bx_new,1,bspace,1)
    call get_time(t1)
    call matvec(n,n_max,space,aspace)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    if (abs(shift).gt.num_thresh) call daxpy(n*n_max,shift,space,1,aspace,1)
    call dgemm('t','n',n_max,n_max,n,one,space,n,aspace,n,zero,a_red,lda)
    call get_time(t1)
    call dsyev('v','l',n_max,a_red,lda,e_red,work,lwork,info)
    call get_time(t2)
    t_diag = t_diag + t2 - t1
    eig = e_red(1:n_max)
!
!   get the ritz vectors:
!
    call dgemm('n','n',n,n_max,n_max,one,space,n,a_red,lda,zero,evec,n)
    call dcopy(n*n_max,evec,1,space,1)
    call dgemm('n','n',n,n_max,n_max,one,aspace,n,a_red,lda,zero,evec,n)
    call dcopy(n*n_max,evec,1,aspace,1)
!
!   if required, also get b times the ritz vector:
!
    if (generalized) then 
      call dgemm('n','n',n,n_max,n_max,one,bspace,n,a_red,lda,zero,evec,n)
      call dcopy(n*n_max,evec,1,bspace,1)
    end if
!
!   do the first iteration explicitly. 
!   build the residuals:
!
    call dcopy(n*n_max,aspace,1,residuals,1)
    if (generalized) then 
      do i_eig = 1, n_max
        call daxpy(n,-eig(i_eig),bspace(:,i_eig),1,residuals(:,i_eig),1)
      end do
    else
      do i_eig = 1, n_max
        call daxpy(n,-eig(i_eig),space(:,i_eig),1,residuals(:,i_eig),1)
      end do
    end if
!
!   compute the preconditioned residuals:
!
    ind_x = 1
    ind_w = ind_x + n_max 
    call precnd(n,n_max,shift-eig(ind_x),residuals(1,ind_x),space(1,ind_w))
!
!   orthogonalize:
!
    call get_time(t1)
    if (generalized) then
      call b_ortho_vs_x(n,n_max,n_max,space,bspace,space(1,ind_w))
!
!     after b_ortho, w is b-orthogonal to x, and orthonormal. 
!     compute the application of b to w, and b-orthonormalize it.
!
      call metvec(n,n_max,space(1,ind_w),bspace(1,ind_w))
      call b_ortho(n,n_max,space(1,ind_w),bspace(1,ind_w))
    else
      call ortho_vs_x(n,n_max,n_max,space,space(1,ind_w),xx,xx)
    end if
    call get_time(t2)
    t_ortho = t_ortho + t2 - t1
!
!   we are now ready to start the main loop.
!   initialize a few parameters
!
    tol_rms = tol
    tol_max = ten*tol
    sqrtn   = sqrt(real(n,dp))
    ok      = .false.
    done    = .false.
    n_act   = n_max
!
    1030 format(t5,'LOBPCG iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    do it = 1, max_iter
!
!     perform the matrix-vector multiplication for this iteration:
!
      call get_time(t1)
      call matvec(n,n_act,space(1,ind_w),aspace(1,ind_w))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
      if (abs(shift).gt.num_thresh) call daxpy(n*n_act,shift,space(1,ind_w),1,aspace(1,ind_w),1)
!
!     build the reduced matrix and diagonalize it:
!
      ld_current = n_max + 2*n_act
      if (it.eq.1) ld_current = 2*n_max
      call dgemm('t','n',ld_current,ld_current,n,one,space,n,aspace,n,zero,a_red,lda)
!
      call get_time(t1)
      call dsyev('v','l',ld_current,a_red,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     if dsyev failed, print an error message and abort (this should not happen)
!
      if (info.ne.0) then
        write(6,'(t3,a,i6)') 'dsyev failed. info = ',info
        stop
      end if
      eig = e_red(1:n_max)
!
!     update x and ax, and, if required, bx:
!
      call dgemm('n','n',n,n_max,ld_current,one,space,n,a_red,lda,zero,x_new,n)
      call dgemm('n','n',n,n_max,ld_current,one,aspace,n,a_red,lda,zero,ax_new,n)
      if (generalized) then
        call dgemm('n','n',n,n_max,ld_current,one,bspace,n,a_red,lda,zero,bx_new,n)
      end if
!
!     compute the residuals and their rms and sup norms:
!
      call dcopy(n*n_max,ax_new,1,residuals,1)
      do i_eig = 1, n_max
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        if (generalized) then
          call daxpy(n,-eig(i_eig),bx_new(:,i_eig),1,residuals(:,i_eig),1)
        else
          call daxpy(n,-eig(i_eig),x_new(:,i_eig),1,residuals(:,i_eig),1)
        end if
        r_norm(1,i_eig) = dnrm2(n,residuals(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(residuals(:,i_eig)))
      end do
!
!     only lock the first converged eigenvalues/vectors.
!
      do i_eig = 1, n_max
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        if (.not. done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information and check for convergence:
!
      if (verbose) then
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*)
      end if
      if (all(done(1:n_targ))) then
        call dcopy(n*n_max,x_new,1,evec,1)
        ok = .true.
        exit
      end if
!
!     compute the number of active eigenvalues. 
!     converged eigenvalues and eigenvectors will be locked and kept 
!     for orthogonalization purposes.
!
      n_act = n_max - count(done)
      ind_x = n_max - n_act + 1
      ind_p = ind_x + n_act
      ind_w = ind_p + n_act
!
!     compute the new p and ap vectors only for the active eigenvectors.
!     this is done by computing the expansion coefficients u_p of x_new
!     -x in the basis of (x,p,w), and then by orthogonalizing then to
!     the coefficients u_x of x_new. 
!
      call mallocate(ld_current,n_max,u_x)
      call mallocate(ld_current,n_act,u_p)
!
      call get_coeffs(lda,ld_current,n_max,n_act,a_red,u_x,u_p)
!
!     p  = space  * u_p
!     ap = aspace * u_p
!     bp = bspace * u_p
!     note that this is numerically safe, as u_p is orthogonal.
!
      call dgemm('n','n',n,n_act,ld_current,one,space,n,u_p,ld_current,zero,evec,n)
      call dcopy(n_act*n,evec,1,space(1,ind_p),1)
      call dgemm('n','n',n,n_act,ld_current,one,aspace,n,u_p,ld_current,zero,evec,n)
      call dcopy(n_act*n,evec,1,aspace(1,ind_p),1)
!
      if (generalized) then
        call dgemm('n','n',n,n_act,ld_current,one,bspace,n,u_p,ld_current,zero,evec,n)
        call dcopy(n_act*n,evec,1,bspace(1,ind_p),1)
      end if
!
      call mfree(u_x)
      call mfree(u_p)
!
!     now, move x_new and ax_new into space and aspace.
!
      call dcopy(n*n_max,x_new,1,space,1)
      call dcopy(n*n_max,ax_new,1,aspace,1)
      if (generalized) then 
        call dcopy(n*n_max,bx_new,1,bspace,1)
      end if
!
!     compute the preconditioned residuals w:
!
      call precnd(n,n_act,shift-eig(1),residuals(1,ind_x),space(1,ind_w))
!
!     orthogonalize w against x and p, and then orthonormalize it:
!
      call get_time(t1)
      if (generalized) then 
        call b_ortho_vs_x(n,n_max+n_act,n_act,space,bspace,space(1,ind_w))
        call metvec(n,n_act,space(1,ind_w),bspace(1,ind_w))
        call b_ortho(n,n_act,space(1,ind_w),bspace(1,ind_w))
      else
        call ortho_vs_x(n,n_max+n_act,n_act,space,space(1,ind_w),xx,xx)
      end if
      call get_time(t2)
      t_ortho = t_ortho + t2 - t1
!
    end do
!
!   deallocate memory and return.
!
    call mfree(work)
    call mfree(tau)
    call mfree(space)
    call mfree(aspace)
    call mfree(residuals)
    call mfree(a_red)
    call mfree(e_red)
    call mfree(x_new)
    call mfree(ax_new)
    call mfree(done)
    call mfree(r_norm)
    if (generalized) then
      call mfree(bspace)
      call mfree(bx_new)
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
    1000 format(t3,'timings for LOBPCG (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
!    

  end subroutine lobpcg_driver
!
  subroutine get_coeffs(lda,ld_current,n_max,n_act,a_red,u_x,u_p)
    implicit none
!
!   given the eigenvetors of the reduced matrix in a_red, extract
!   the expansion coefficients for x_new (u_x) and assemble the 
!   ones for p_new in u_p.
!
!   the coefficients u_p are computed as the difference between the
!   coefficients for x_new and x_old, and only the columns associated
!   with active eigenvectors are considered. 
!   u_p is then orthogonalized to u_x: this not only guarantees that
!   the p_new vectors will be orthogonal to x_new, but also allows one
!   to reuse the ax, aw, and ap vectors to compute ap_new, without
!   loosing numerical precision.
!
    integer,                          intent(in)    :: lda, ld_current, n_max, n_act
    real(dp), dimension(lda,lda), intent(in)    :: a_red
    real(dp), dimension(ld_current,n_max), intent(inout) :: u_x
    real(dp), dimension(ld_current,n_act), intent(inout) :: u_p
!
    integer               :: ind_x, off_x, i_eig
    real(dp)              :: xx(1)
!
    off_x = n_max - n_act
    ind_x = off_x + 1
!
    u_x(1:ld_current,1:n_max) = a_red(1:ld_current,1:n_max)
!
!   u_p = u_x for the active vectors only
!
    u_p = u_x(:,ind_x:n_max)
!
!   remove the coefficients for x from u_p
!
    do i_eig = 1, n_act
      u_p(off_x + i_eig,i_eig) = u_p(off_x + i_eig,i_eig) - one
    end do
!
!   orthogonalize:
!
    call ortho_vs_x(ld_current,n_max,n_act,u_x,u_p,xx,xx)
!
!   all done.
!
    return
!
  end subroutine get_coeffs

end module mod_lobpcg_driver

