module diaglib
  use real_precision
  implicit none
!
! diaglib - a fortran library of matrix-free iterative algorithms to
! compute a few eigenvalues and eigenvectors of large matrices.
! ==================================================================
!
! Implementation by
!
!   Ivan Gianni', Tommaso Nottoli, and Filippo Lipparini
!   MoLECoLab Pisa
!   Department of Chemistry and Industrial Chemistry
!   University of Pisa
!   Via G. Moruzzi 13, I-56124, Pisa, Italy
!
! Pisa, november 2022
!
!                                                                     
!                                       mm                            
!                                    mMMm                             
!                                  mMMMMm         m                   
!                                 mMMMMm          mMm                 
!                                 mMMMMm          mMm                 
!                                 mMMMMMm        mMMm                 
!                                 MMMMMMMMMMMMMMMMMMm                 
!                                mMMMMMMMMMMMMMMMMMm                  
!       __  ___      __    ____________      __MMMm     __            
!      /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_           
!     / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \          
!    / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /          
!   /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/           
!           /_  __/ __ \/ __ \/ /  / ___/                             
!            / / / / / / / / / /   \__ \                              
!           / / / /_/ / /_/ / /___ __/ /                              
!          /_/  \____/\____/_____/____/                               
!            mMMMMMMMMMMMMMMMm                                        
!          mMMMMMMMMMMMMMMMm                                          
!        mMMMMMMMMMMMMMMMMM   + ------------------------------------ +
!       mMMMMMMMMMMMMMMMMm    |            D I A G L I B             |
!      mMMMMMMMMMMMMMMMMMm    + ------------------------------------ +
!      mMMMMm       mMMMMMm   | I. Gianni', T. Nottoli, F. Lipparini |
!      mMMMm       mMMMMMMm   |                                      |
!       mMm       mMMMMMMm    |                            ver 1.0   |
!        m       mMMMMMMm     |              molecolab.dcci.unipi.it |
!               mMMMMMm       + ------------------------------------ +
!                                                                     
! description of the library:
! ===========================
!
! diaglib provides an implementation of two matrix-free algorithms to
! compute a few eigenvalues and eigenvectors of a large, possibly sparse
! matrix.
!
! the available algorithms are
!
! 1) locally optimal block preconditioned conjugate gradient 
!
! 2) davidson-liu
!
! both algorithms require two user-provided routines to apply the matrx
! and a suitable preconditioner to a set of vectors.
! such routines have the following interface:
!
!   subroutine matvec(n,m,x,ax)
!   subroutine precnd(n,m,shift,x,ax)
!
! where n,m are integers and x(n,m) and ax(n,m) are double precision
! arrays.
! as using the first eigenvalue in a shift-and-invert spirit is very 
! common, a double precision scalar shift is also passed to precnd.
!
! both implementations favor numerical stability over efficiency and are
! targeted at applications in molecular quantum chemistry, such as in
! (full) ci or augmented hessian calculations, where typically m << n.
!
! list of provided routines:
! ==========================
!
! lobpcg_driver:   main lobpcg driver.
!
! davidson_driver: main davidson-liu driver
!
! ortho_vs_x:      subroutine to orthogonalize a set of vectors w against
!                  another set of orthonormal vectors v.
!                  the resulting w vectors are then also orthonormalized.
!                  the linear transformation applied to w can also be 
!                  applied to a second set of vectors, tipically, aw.
!
! ortho:           subroutine to perferm the orthonormalization of a set
!                  of vectors v using QR decomposition.
!                  the linear transformation applied to v can also be 
!                  applied to a second set of vectors, tipically, av.
!
! ortho_cd:        subroutine to perferm the orthonormalization of a set
!                  of vectors v using cholesky decomposition with level
!                  shifting and iterative refinement.
!                  the linear transformation applied to v can also be 
!                  applied to a second set of vectors, tipically, av.
!
! b_ortho_vs_x:    same as ortho_vs_x, but vectors are b-orthogonalized
!                  against the existing ones, where b is the metric in the
!                  generalized eigenvalues problem.
!
! b_ortho:         same as ortho_cd, but the orthogonalization is performed
!                  with respect to the metric b.
!
! check_guess:     routine to check whether a guess is provided by the user
!                  and whether the provided vectors are already orthonormal.
!                  if no guess is provided, a random guess is generated.
!
! get_coeffs:      utility to compute the expasion coefficient for the p
!                  vectors. used in lobpcg.
!
! diag_shift:      utility to apply a diagonal shift to the diagonal of a
!                  given matrix. 
!                  used in ortho_cd.
!
! get_mem_lapack:  utility to compute how much scratch memory is required by
!                  the various lapack routines called in the code.
!
! check_mem:       utility to check whether memory was allocated or deallocated
!                  successfully.
!
! prtmat:          utility for debug printing.
!
! get_time:        subroutine to get cpu and wall time.
!
! please, see the specific routines for a more detailed description.
!
! dependencies:
! =============
!
! the following blas/lapack routines are used:
!
!   blas:          daxpy, dcopy, dgemm, dnrm2
!   lapack:        dgeqrf, dpotrf, dsyev, dtrsm 
!
! shared variables:
! =================
!
! useful constants
!
  real(dp), parameter    :: zero = 0.0_dp, one = 1.0_dp, ten = 10.0_dp
!
! convergence thresholds for orthogonalization
!
  real(dp), parameter    :: tol_ortho = 1.0e-14_dp
! 
! memory and info for lapack routines
!
  integer                :: lwork, info
  real(dp), allocatable  :: work(:), tau(:)
!
! timings:
!
  real(dp)               :: t1(2), t2(2), t_diag(2), t_ortho(2), &
                            t_mv(2), t_tot(2)
!
! subroutines:
! ============
!
  contains
!
  subroutine lobpcg_driver(verbose,gen_eig,n,n_targ,n_max,max_iter,tol, &
                           shift,matvec,precnd,bvec,eig,evec,ok)
!
!   main driver for lobpcg.
!
!   input variables:
!   ================
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   gen_eig:  logical, whether a generalized eigenvalue problem has
!             to be solved. 
!
!   n:        integer, size of the matrix to be diagonalized.
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   shift:    double precision real, a diagonal level shifting parameter
!
!   matvec:   external subroutine that performs the matrix-vector
!             multiplication
!
!   precnd:   external subroutine that applies a preconditioner.
!
!   bvec:     external subroutine that applies the metric to a vector.
!             only referenced is gen is true.
!
!   output variables:
!   =================
!
!   eig:      double precision array of size n_max. if ok is true, 
!             the computed eigenvalues in asceding order.
!
!   evec:     double precision array of size (n,n_max). 
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if lobpcg converged.
!
    logical,                      intent(in)    :: verbose, gen_eig
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                      intent(inout) :: ok
    external                                    :: matvec, precnd, bvec
!
!   local variables:
!   ================
!
    integer               :: it, i_eig, n_act, ind_x, ind_w, ind_p, &
                             len_a, len_u
    integer               :: istat
    real(dp)              :: sqrtn, tol_rms, tol_max, xx(1)
    real(dp), allocatable :: space(:,:), aspace(:,:), bspace(:,:), &
                             a_red(:,:), e_red(:), r(:,:), r_norm(:,:)
    real(dp), allocatable :: u_x(:,:), u_p(:,:), x_new(:,:), ax_new(:,:), &
                             bx_new(:,:)
    logical,  allocatable :: done(:)
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dnrm2, daxpy, dsyev, dcopy
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(2*n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    len_a = 3*n_max
    allocate (space(n,len_a), aspace(n,len_a), bspace(n,len_a), &
              r(n,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(len_a,len_a), e_red(len_a), stat=istat)
    call check_mem(istat)
!
!   allocate memory for temporary copies of x, ax, and bx:
!
    allocate (x_new(n,n_max), ax_new(n,n_max), bx_new(n,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   clean out:
!
    t_diag   = zero
    t_ortho  = zero
    t_mv     = zero
    t_tot    = zero
    space    = zero
    bspace   = zero
    aspace   = zero
    a_red    = zero
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
    if (gen_eig) then 
      call bvec(n,n_max,evec,bx_new)
      call b_ortho(n,n_max,evec,bx_new)
    end if
!
!   compute the first eigenpairs by diagonalizing the reduced matrix:
!
    call dcopy(n*n_max,evec,1,space,1)
    if (gen_eig) call dcopy(n*n_max,bx_new,1,bspace,1)
    call get_time(t1)
    call matvec(n,n_max,space,aspace)
    call get_time(t2)
    t_mv = t_mv + t2 - t1
    if (shift.ne.zero) call daxpy(n*n_max,shift,space,1,aspace,1)
    call dgemm('t','n',n_max,n_max,n,one,space,n,aspace,n,zero,a_red,len_a)
    call get_time(t1)
    call dsyev('v','l',n_max,a_red,len_a,e_red,work,lwork,info)
    call get_time(t2)
    t_diag = t_diag + t2 - t1
    eig = e_red(1:n_max)
!
!   get the ritz vectors:
!
    call dgemm('n','n',n,n_max,n_max,one,space,n,a_red,len_a,zero,evec,n)
    call dcopy(n*n_max,evec,1,space,1)
    call dgemm('n','n',n,n_max,n_max,one,aspace,n,a_red,len_a,zero,evec,n)
    call dcopy(n*n_max,evec,1,aspace,1)
!
!   if required, also get b times the ritz vector:
!
    if (gen_eig) then 
      call dgemm('n','n',n,n_max,n_max,one,bspace,n,a_red,len_a,zero,evec,n)
      call dcopy(n*n_max,evec,1,bspace,1)
    end if
!
!   do the first iteration explicitly. 
!   build the residuals:
!
    call dcopy(n*n_max,aspace,1,r,1)
    if (gen_eig) then 
      do i_eig = 1, n_max
        call daxpy(n,-eig(i_eig),bspace(:,i_eig),1,r(:,i_eig),1)
      end do
    else
      do i_eig = 1, n_max
        call daxpy(n,-eig(i_eig),space(:,i_eig),1,r(:,i_eig),1)
      end do
    end if
!
!   compute the preconditioned residuals:
!
    ind_x = 1
    ind_w = ind_x + n_max 
    call precnd(n,n_max,shift-eig(ind_x),r(1,ind_x),space(1,ind_w))
!
!   orthogonalize:
!
    call get_time(t1)
    if (gen_eig) then
      call b_ortho_vs_x(n,n_max,n_max,space,bspace,space(1,ind_w))
!
!     after b_ortho, w is b-orthogonal to x, and orthonormal. 
!     compute the application of b to w, and b-orthonormalize it.
!
      call bvec(n,n_max,space(1,ind_w),bspace(1,ind_w))
      call b_ortho(n,n_max,space(1,ind_w),bspace(1,ind_w))
    else
      call ortho_vs_x(.false.,n,n_max,n_max,space,space(1,ind_w),xx,xx)
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
      if (shift.ne.zero) call daxpy(n*n_act,shift,space(1,ind_w),1,aspace(1,ind_w),1)
!
!     build the reduced matrix and diagonalize it:
!
      len_u = n_max + 2*n_act
      if (it.eq.1) len_u = 2*n_max
      call dgemm('t','n',len_u,len_u,n,one,space,n,aspace,n,zero,a_red,len_a)
!
      call get_time(t1)
      call dsyev('v','l',len_u,a_red,len_a,e_red,work,lwork,info)
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
      call dgemm('n','n',n,n_max,len_u,one,space,n,a_red,len_a,zero,x_new,n)
      call dgemm('n','n',n,n_max,len_u,one,aspace,n,a_red,len_a,zero,ax_new,n)
      if (gen_eig) then
        call dgemm('n','n',n,n_max,len_u,one,bspace,n,a_red,len_a,zero,bx_new,n)
      end if
!
!     compute the residuals and their rms and sup norms:
!
      call dcopy(n*n_max,ax_new,1,r,1)
      do i_eig = 1, n_max
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        if (gen_eig) then
          call daxpy(n,-eig(i_eig),bx_new(:,i_eig),1,r(:,i_eig),1)
        else
          call daxpy(n,-eig(i_eig),x_new(:,i_eig),1,r(:,i_eig),1)
        end if
        r_norm(1,i_eig) = dnrm2(n,r(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(r(:,i_eig)))
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
        call dcopy(n*n_max,ax_new,1,evec,1)
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
      allocate (u_x(len_u,n_max), u_p(len_u,n_act), stat = istat)
      call check_mem(istat)
!
      call get_coeffs(len_a,len_u,n_max,n_act,a_red,u_x,u_p)
!
!     p  = space  * u_p
!     ap = aspace * u_p
!     bp = bspace * u_p
!     note that this is numerically safe, as u_p is orthogonal.
!
      call dgemm('n','n',n,n_act,len_u,one,space,n,u_p,len_u,zero,evec,n)
      call dcopy(n_act*n,evec,1,space(1,ind_p),1)
      call dgemm('n','n',n,n_act,len_u,one,aspace,n,u_p,len_u,zero,evec,n)
      call dcopy(n_act*n,evec,1,aspace(1,ind_p),1)
!
      if (gen_eig) then
        call dgemm('n','n',n,n_act,len_u,one,bspace,n,u_p,len_u,zero,evec,n)
        call dcopy(n_act*n,evec,1,bspace(1,ind_p),1)
      end if
!
      deallocate(u_x, u_p, stat = istat)
      call check_mem(istat)
!
!     now, move x_new and ax_new into space and aspace.
!
      call dcopy(n*n_max,x_new,1,space,1)
      call dcopy(n*n_max,ax_new,1,aspace,1)
      if (gen_eig) then 
        call dcopy(n*n_max,bx_new,1,bspace,1)
      end if
!
!     compute the preconditioned residuals w:
!
      call precnd(n,n_act,shift-eig(1),r(1,ind_x),space(1,ind_w))
!
!     orthogonalize w against x and p, and then orthonormalize it:
!
      call get_time(t1)
      if (gen_eig) then 
        call b_ortho_vs_x(n,n_max+n_act,n_act,space,bspace,space(1,ind_w))
        call bvec(n,n_act,space(1,ind_w),bspace(1,ind_w))
        call b_ortho(n,n_act,space(1,ind_w),bspace(1,ind_w))
      else
        call ortho_vs_x(.false.,n,n_max+n_act,n_act,space,space(1,ind_w),xx,xx)
      end if
      call get_time(t2)
      t_ortho = t_ortho + t2 - t1
!
    end do
!
    call get_time(t2)
    t_tot = t2 - t_tot
!
!   if required, print timings
!
    1000 format(t3,'timings for lobpcg (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!
!   deallocate memory and return.
!
    deallocate (work, tau, space, aspace, bspace, r, a_red, e_red, & 
                x_new, ax_new, bx_new, done, r_norm, stat = istat)
    call check_mem(istat)
!
    return

  end subroutine lobpcg_driver
!
  subroutine davidson_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav, &
                             shift,matvec,precnd,eig,evec,ok)
!
!   main driver for davidson-liu.
!
!   input variables:
!   ================
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   n:        integer, size of the matrix to be diagonalized.
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   max_dav:  integer. maximum allowed number of iterations before a 
!             restart is forced. 
!             when n_max eigenvalues are searched, this implies a maximum
!             dimension of the expansion subspace of n_max * max_dav.
!
!   shift:    double precision real, a diagonal level shifting parameter
!
!   matvec:   external subroutine that performs the matrix-vector
!             multiplication
!
!   precnd:   external subroutine that applied a preconditioner.
!
!   output variables:
!   =================
!
!   eig:      double precision array of size n_max. if ok is true, 
!             the computed eigenvalues in asceding order.
!
!   evec:     double precision array of size (n,n_max). 
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if davidson converged.
!
    logical,                      intent(in)    :: verbose
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter, max_dav
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                      intent(inout) :: ok
    external                                    :: matvec, precnd
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: it, i_eig
!
    real(dp)              :: sqrtn, tol_rms, tol_max
    real(dp)              :: xx(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: space(:,:), aspace(:,:), r(:,:), r_norm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), e_red(:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dsyev
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (space(n,lda), aspace(n,lda), r(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda,lda), a_copy(lda,lda), e_red(lda), stat=istat)
    call check_mem(istat)
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
!   initialize the number of active vectors and the associated indices.
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
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
      ldu = ldu + n_act
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
      call dgemm('t','n',ldu,n_act,n,one,space,n,aspace(1,i_beg+n_rst),n,zero,a_red(1,i_beg+n_rst),lda)
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
      call dsyev('v','u',ldu,a_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      eig = e_red(1:n_max)
!
      call dgemm('n','n',n,n_max,ldu,one,space,n,a_copy,lda,zero,evec,n)
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,aspace,n,a_copy,lda,zero,r,n)
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),evec(:,i_eig),1,r(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,r(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(r(:,i_eig)))
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
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        m_dim = m_dim + 1
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
        call precnd(n,n_act,-eig(ind),r(1,ind),space(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call ortho_vs_x(.false.,n,ldu,n_act,space,space(1,i_beg),xx,xx)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        n_act = n_max
        space = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        call dcopy(n_max*n,evec,1,space,1)
        aspace = zero
        a_red  = zero
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1 
        m_dim = 1
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
    call get_time(t2)
    t_tot = t2 - t_tot
!
!   if required, print timings
!
    1000 format(t3,'timings for davidson (cpu/wall): ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!      
!   deallocate memory
!
    deallocate (work, tau, space, aspace, r, done, r_norm, a_red, a_copy, e_red)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine davidson_driver
!
  subroutine gen_david_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav, &
                              shift,matvec,precnd,bvec,eig,evec,ok)
!
!   main driver for davidson-liu.
!
!   input variables:
!   ================
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   n:        integer, size of the matrix to be diagonalized.
!
!   n_targ:   integer, number of required eigenpairs.
!
!   n_max:    integer, maximum size of the search space. should be 
!             >= n_targ. note that eig and evec should be allocated 
!             n_max and (n,n_max) rather than n_targ and (n,n_targ). 
!             for better convergence, a value larger than n_targ (eg.,
!             n_targ + 10) is recommended.
!   
!   max_iter: integer, maximum allowed number of iterations.
!
!   tol:      double precision real, the convergence threshold.
!
!   max_dav:  integer. maximum allowed number of iterations before a 
!             restart is forced. 
!             when n_max eigenvalues are searched, this implies a maximum
!             dimension of the expansion subspace of n_max * max_dav.
!
!   shift:    double precision real, a diagonal level shifting parameter
!
!   matvec:   external subroutine that performs the matrix-vector
!             multiplication
!
!   precnd:   external subroutine that applies a preconditioner.
!
!   bvec:     external subroutine that applies the metric.
!
!   output variables:
!   =================
!
!   eig:      double precision array of size n_max. if ok is true, 
!             the computed eigenvalues in asceding order.
!
!   evec:     double precision array of size (n,n_max). 
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if davidson converged.
!
    logical,                      intent(in)    :: verbose
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter, max_dav
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                      intent(inout) :: ok
    external                                    :: matvec, precnd, bvec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: m_dim, ldu
!
!   number of frozen (i.e. converged) vectors
!
    integer               :: n_frozen
!
    integer               :: it, i_eig
!
    real(dp)              :: sqrtn, tol_rms, tol_max
    real(dp)              :: xx(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: space(:,:), aspace(:,:), bspace(:,:), r(:,:), r_norm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
!
!   scratch:
!
    real(dp), allocatable :: b_evec(:,:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dcopy, dnrm2, dgemm, dsyev
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (space(n,lda), aspace(n,lda), bspace(n,lda), r(n,n_max), &
              b_evec(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda,lda), a_copy(lda,lda), e_red(lda), s_red(lda,lda), &
              s_copy(lda,lda), stat=istat)
    call check_mem(istat)
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
    bspace  = zero
    a_red   = zero
    s_red   = zero
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
!   initialize the number of active vectors and the associated indices.
!
    n_act = n_max
    ind   = 1
    i_beg = 1
!
!   initialize the counter for the expansion of the subspace
!
    m_dim = 1
    ldu   = 0
!
!   initialize to false the restart
!
    restart = .false.
!
!   main loop:
!
    1030 format(t5,'Generalized Davidson-Liu iterations (tol=',d10.2,'):',/, &
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
      ldu = ldu + n_act
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call matvec(n,n_act,space(1,i_beg+n_rst),aspace(1,i_beg+n_rst))
      call bvec(n,n_act,space(1,i_beg+n_rst),bspace(1,i_beg+n_rst))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrices 
!
      call dgemm('t','n',ldu,n_act,n,one,space,n,aspace(1,i_beg+n_rst),n,zero,a_red(1,i_beg+n_rst),lda)
      call dgemm('t','n',ldu,n_act,n,one,space,n,bspace(1,i_beg+n_rst),n,zero,s_red(1,i_beg+n_rst),lda)
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
      s_copy = s_red
!
!     diagonalize the reduced matrix
!
      call get_time(t1)
      call dsygv(1,'v','u',ldu,a_copy,lda,s_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      eig = e_red(1:n_max)
!
      call dgemm('n','n',n,n_max,ldu,one,space,n,a_copy,lda,zero,evec,n)
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,aspace,n,a_copy,lda,zero,r,n)
      call dgemm('n','n',n,n_max,ldu,one,bspace,n,a_copy,lda,zero,b_evec,n)
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),b_evec(:,i_eig),1,r(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,r(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(r(:,i_eig)))
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
      if (m_dim .lt. dim_dav) then
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
        m_dim = m_dim + 1
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
        call precnd(n,n_act,-eig(ind),r(1,ind),space(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call ortho_vs_x(.false.,n,ldu,n_act,space,space(1,i_beg),xx,xx)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        n_act = n_max
        space = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        call dcopy(n_max*n,evec,1,space,1)
        call dcopy(n_max*n,b_evec,1,bspace,1)
        aspace = zero
        bspace = zero
        a_red  = zero
        s_red  = zero
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1 
        m_dim = 1
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
    call get_time(t2)
    t_tot = t2 - t_tot
!
!   if required, print timings
!
    1000 format(t3,'timings for davidson (cpu/wall): ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!      
!   deallocate memory
!
    deallocate (work, tau, space, aspace, r, done, r_norm, a_red, a_copy, e_red)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine gen_david_driver
!
! orthogonalization routines:
! ===========================
!
  subroutine ortho(do_other,n,m,u,w)
    implicit none
!
!   orthogonalize m vectors of lenght n using the QR decomposition.
!   this is done by computing U = QR and then by solving the upper 
!   triangular system U(ortho)R = U.
!
!   using this strategy allows to apply the same linear transformation
!   that orthogonalizes U to a second set of vectors that usually contain
!   the product AU, where A is some matrix that has already been applied
!   to U. 
!   this is useful when U and AU are built together without explicitly 
!   performing the matrix vector multiplication.
!
!   arguments:
!   ==========
!
    logical,                     intent(in)    :: do_other
    integer,                     intent(in)    :: n, m
    real(dp),  dimension(n,m),   intent(inout) :: u, w
!
!   local scratch
!   =============
!
    real(dp), allocatable :: v(:,:)
!
!   external functions:
!   ===================
!
    external dgeqrf, dtrsm
!
    allocate (v(n,m))
    v = u
    call dgeqrf(n,m,u,n,tau,work,lwork,info)
!
    call dtrsm('r','u','n','n',n,m,one,u,n,v,n)
!
    if (do_other) call dtrsm('r','u','n','n',n,m,one,u,n,w,n)
!
    u = v
!
    deallocate (v)
    return
  end subroutine ortho
!
  subroutine b_ortho(n,m,u,bu)
    implicit none
!
!   b-orthogonalize m vectors of lenght n using the cholesky decomposition
!   of the overlap matrix.
!   this is in principle not a good idea, as the u'bu matrix can be very
!   ill-conditioned, independent of how bas is b, and only works if x is
!   already orthonormal. 
!
!   arguments:
!   ==========
!
    integer,                     intent(in)    :: n, m
    real(dp),  dimension(n,m),   intent(inout) :: u, bu
!
!   local variables
!   ===============
!
    integer               :: info, i, j
    real(dp), allocatable :: metric(:,:), sigma(:), u_svd(:,:), vt_svd(:,:), &
                             temp(:,:)
    real(dp), parameter   :: tol_svd = 1.0e-5_dp
    logical,  parameter   :: use_svd = .true.
!
!   external functions:
!   ===================
!
    external dpotrf, dtrsm, dgemm
!
    allocate (metric(m,m))
!
    call dgemm('t','n',m,m,n,one,u,n,bu,n,zero,metric,m)
!
    if (use_svd) then
!
!     debug option: use svd to b-orthonormalize, by computing
!     b**(-1/2)
!
      allocate (sigma(m), u_svd(m,m), vt_svd(m,m), temp(n,m))
      call dgesvd('a','a',m,m,metric,m,sigma,u_svd,m,vt_svd,m,work,lwork,info)
!
!     compute sigma**(-1/2)
!
      do i = 1, m
        if (sigma(i) .gt. tol_svd) then
          sigma(i) = 1/sqrt(sigma(i))
        else
          sigma(i) = zero
        end if
      end do
!
!     compute metric ** (-1/2). first, compute sigma ** (-1/2) vt
!
      metric = zero
      do i = 1, m
        do j = 1, m
          metric(j,i) = metric(j,i) + sigma(j)*vt_svd(j,i)
        end do
      end do
!
!     now, multiply for u:
!
      vt_svd = metric
      call dgemm('n','n',m,m,m,one,u_svd,m,vt_svd,m,zero,metric,m)
!
!     metric contains s ** (-1/2), and projects out directions corresponding
!     to pathological singular values. 
!     orthogonalize u and bu:
!
      call dgemm('n','n',n,m,m,one,u,n,metric,m,zero,temp,n)
      u = temp
      call dgemm('n','n',n,m,m,one,bu,n,metric,m,zero,temp,n)
      bu = temp
!
      deallocate (sigma, u_svd, vt_svd, temp)
    else
!
!     compute the cholesky factorization of the metric.
!
      call dpotrf('l',m,metric,m,info)
!
!     get u * l^-T and bu * l^-T
!
      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,bu,n)
    end if
!
    deallocate (metric)
    return
  end subroutine b_ortho
!
  subroutine ortho_cd(do_other,n,m,u,w,ok)
    implicit none
!
!   orthogonalize m vectors of lenght n using the Cholesky factorization
!   of their overlap. 
!   this is done by metric = U^t U and then by computing its cholesky 
!   decompositoin metric = L L^t. The orthogonal vectors are obtained then
!   by solving the triangular linear system
!
!     U(ortho)L^T = U
!
!   using this strategy allows to apply the same linear transformation
!   that orthogonalizes U to a second set of vectors that usually contain
!   the product AU, where A is some matrix that has already been applied
!   to U.  this is useful when U and AU are built together without explicitly 
!   performing the matrix vector multiplication.
!
!   as cholesky decomposition is not the most stable way of orthogonalizing
!   a set of vectors, the orthogonalization is refined iteratively. 
!   still, this routine can fail. 
!   a logical flag is then set to false, so that the calling program can 
!   call a more robust orthogonalization routine without aborting.
!
!   arguments:
!   ==========
!
    logical,                   intent(in)    :: do_other
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(inout) :: u, w
    logical,                   intent(inout) :: ok
!
!   local variables
!   ===============
!
    integer               :: it, it_micro
    real(dp)              :: metric_norm, dnrm2, alpha, unorm, shift
    logical               :: macro_done, micro_done
    integer, parameter    :: maxit = 10
!
!   local scratch
!   =============
!
    real(dp), allocatable :: metric(:,:), msave(:,:)
!
!   external functions:
!   ===================
!
    external              :: dgemm, dgeqrf, dtrsm, dnrm2
!
!   get memory for the metric.
!
    allocate (metric(m,m), msave(m,m))
    macro_done = .false.
!
!   assemble the metric
!
    it = 0
    do while(.not. macro_done)
      it = it + 1
      if (it .gt. maxit) then
!
!       ortho_cd failed. return with an error message
!
        ok = .false.
        write(6,100) ' maximum number of iterations reached.'
        return
      end if
      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
!
!   compute the cholesky factorization of the metric.
!
      call dpotrf('l',m,metric,m,info)
!
!     if dpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!
      if (info.ne.0) then
!
        alpha      = 100.0_dp
        unorm      = dnrm2(n*m,u,1)
        it_micro   = 0
        micro_done = .false.
!
!       add larger and larger shifts to the diagonal until dpotrf manages to factorize it.
!
        do while (.not. micro_done)
          it_micro = it_micro + 1
          if (it_micro.gt.maxit) then
!
!           something went very wrong. return with an error status, the orthogonalization
!           will be carried out using a different algorithm.
!
            ok = .false.
            write(6,100) ' maximum number of iterations for factorization reached.'
            return
          end if
!
          shift = max(epsilon(one)*alpha*unorm,tol_ortho)
          metric = msave
          call diag_shift(m,shift,metric)
          call dpotrf('l',m,metric,m,info)
          alpha = alpha * 10.0_dp
          micro_done = info.eq.0
        end do
!
      end if
!
!     orthogonalize u by computing a solution to u(ortho) l^t = u
!     if required, apply the same transformation to w.
!    
      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
      if (do_other) call dtrsm('r','l','t','n',n,m,one,metric,m,w,n)
!
!     check that the vectors in v are really orthogonal
!
      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
      metric_norm = abs(dnrm2(m*m,metric,1) - sqrt(dble(m)))
      macro_done = metric_norm .lt. tol_ortho
    end do
!
    100 format(t3,'ortho_cd failed with the following error:',a)
!
    ok = .true.
!
    deallocate (metric)
    return
  end subroutine ortho_cd
!
  subroutine ortho_vs_x(do_other,n,m,k,x,u,ax,au)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x 
!   is assumed to be orthogonal, orthogonalize u against x.
!
!   if required, orthogonalize au to ax using the same linear
!   transformation, where ax and au are the results of the 
!   application of a matrix a to both x and u.
!
!   furthermore, orthonormalize u and, if required, apply the
!   same transformation to au.
!
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold. 
!
!   arguments:
!   ==========
!
    logical,                   intent(in)    :: do_other
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x, ax
    real(dp),  dimension(n,k), intent(inout) :: u, au
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it
    real(dp)               :: xu_norm
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   ===================
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
!   
    integer, parameter     :: maxit = 10
    logical, parameter     :: useqr = .true.
!
!   allocate space for the overlap between x and u.
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it   = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
    if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
    if (.not. ok .or. useqr) call ortho(do_other,n,k,u,au)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^t u)
!
      call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
      if (do_other) call dgemm('n','n',n,k,m,-one,ax,n,xu,m,one,au,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
      if (.not. ok .or. useqr) call ortho(do_other,n,k,u,au)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
      xu_norm = dnrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine ortho_vs_x
!
  subroutine b_ortho_vs_x(n,m,k,x,bx,u)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x 
!   is assumed to be orthogonal, b-orthogonalize u against x.
!   furthermore, orthonormalize u.
!
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold. 
!
!   arguments:
!   ==========
!
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x, bx
    real(dp),  dimension(n,k), intent(inout) :: u
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it
    real(dp)               :: xu_norm, xx(1)
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   ===================
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
!   
    integer, parameter     :: maxit = 10
    logical, parameter     :: useqr = .true.
!
!   allocate space for the overlap between x and u.
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it   = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
    if (.not. useqr) call ortho_cd(.false.,n,k,u,xx,ok)
    if (.not. ok .or. useqr) call ortho(.false.,n,k,u,xx)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (bx^t u)
!
      call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd(.false.,n,k,u,xx,ok)
      if (.not. ok .or. useqr) call ortho(.false.,n,k,u,xx)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
      xu_norm = dnrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of b_ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine b_ortho_vs_x
!
! utilities
! =========
!
  subroutine diag_shift(n,shift,a)
    implicit none
!  
!   add shift to the diagonal elements of the matric a
!  
    integer,                  intent(in)    :: n
    real(dp),                 intent(in)    :: shift
    real(dp), dimension(n,n), intent(inout) :: a
!  
    integer :: i
!  
    do i = 1, n
      a(i,i) = a(i,i) + shift
    end do
!  
    return
  end subroutine diag_shift
!
  subroutine get_coeffs(len_a,len_u,n_max,n_act,a_red,u_x,u_p)
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
    integer,                          intent(in)    :: len_a, len_u, n_max, n_act
    real(dp), dimension(len_a,len_a), intent(in)    :: a_red
    real(dp), dimension(len_u,n_max), intent(inout) :: u_x
    real(dp), dimension(len_u,n_act), intent(inout) :: u_p
!
    integer               :: ind_x, off_x, i_eig
    real(dp)              :: xx(1)
!
    off_x = n_max - n_act
    ind_x = off_x + 1
!
    u_x(1:len_u,1:n_max) = a_red(1:len_u,1:n_max)
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
    call ortho_vs_x(.false.,len_u,n_max,n_act,u_x,u_p,xx,xx)
!
!   all done.
!
    return
!
  end subroutine get_coeffs
!
  subroutine check_guess(n,m,evec)
    implicit none
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(inout) :: evec
!
    integer               :: i, j, istat
    real(dp)              :: fac, diag_norm, out_norm, xx(1)
!
    real(dp), allocatable :: overlap(:,:)
    real(dp)              :: dnrm2
    external              :: dnrm2
!
!   check whether evec is zero.
!
    fac = dnrm2(n*m,evec,1)
    if (fac.eq.zero) then
!
!     no luck. make a random guess, then orthonormalize it.
!
      call random_number(evec)
      call ortho(.false.,n,m,evec,xx)
    else
!
!     compute the overlap and check that the vectors are orthonormal.
!
      allocate (overlap(m,m), stat=istat)
      call check_mem(istat)
      call dgemm('t','n',m,m,n,one,evec,n,evec,n,zero,overlap,m)
      diag_norm = zero
      out_norm  = zero
      do i = 1, m
        diag_norm = diag_norm + overlap(i,i)**2
        do j = 1, i-1
          out_norm = out_norm + overlap(j,i)**2
        end do
      end do
!
      diag_norm = diag_norm/real(m,dp)
!
      if (diag_norm .ne. one .or. out_norm.ne.zero) then
!
!       orthogonalize the guess:
!
        call ortho(.false.,n,m,evec,xx)
      end if
!
      deallocate (overlap, stat = istat)
      call check_mem(istat)
    end if
!
    return
  end subroutine check_guess
!
!
  subroutine check_mem(istat)
    implicit none
!
!   silly routine that checks that istat is zero after memory allocation
!   and aborts otherwise printing an error message.
!
    integer, intent(in) :: istat
!
    1000 format(t3,'memory allocation failed. istat = ',i8)
    if (istat.ne.0) then
      write(6,1000) istat
      stop
    end if
    return
  end subroutine check_mem
!
  integer function get_mem_lapack(n,n_max)
    integer, intent(in)    :: n, n_max
!
    integer           :: lwork1, lwork2, len_rr, len_qr, nb
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
    get_mem_lapack = max(lwork1,lwork2)
    return
  end function get_mem_lapack
! 
  subroutine prtmat(n,lda,mat,iform)
    integer,                      intent(in) :: n, lda, iform
    real(dp), dimension(lda,lda), intent(in) :: mat
!
    integer :: i
!
    100 format(t3,20f12.6)
    200 format(t3,20d12.4)
!
    do i = 1, n
      if (iform.eq.1) write(6,100) mat(i,:)
      if (iform.eq.2) write(6,200) mat(i,:)
    end do
    return
  end subroutine prtmat
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
end module diaglib
