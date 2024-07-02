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
!   Ivan Gianni', Tommaso Nottoli, Federica Pes, Antoine Levitt, and Filippo Lipparini
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
!      mMMMMm       mMMMMMm   | I. Gianni', T. Nottoli, F. Pes       |
!      mMMMm       mMMMMMMm   | A. Levitt, F. Lipparini              !
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
  private
!
! useful constants
!
  real(dp), parameter    :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
!
! convergence thresholds for orthogonalization
!
  real(dp), parameter    :: tol_ortho = 1.0d-12 ! two * epsilon(one)
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
  public :: lobpcg_driver, davidson_driver, gen_david_driver, caslr_driver, caslr_eff_driver, &
            nonsym_driver, ortho, b_ortho, ortho_cd, ortho_vs_x, b_ortho_vs_x   
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
        call ortho_vs_x(n,n_max+n_act,n_act,space,space(1,ind_w),xx,xx)
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
  subroutine caslr_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                          apbmul,ambmul,spdmul,smdmul,lrprec,eig,evec,ok)
    use utils
    implicit none
    logical, intent(in)                          :: verbose
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n_max),    intent(inout) :: eig
    real(dp), dimension(n2,n_max), intent(inout) :: evec
    logical,                       intent(inout) :: ok
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda, lda2
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
    integer               :: i, j, k, it, i_eig
!
    integer               :: lwork_svd
!
    real(dp)              :: sqrtn, tol_rms, tol_max, growth
    real(dp)              :: xx(1), lw_svd(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    real(dp), allocatable :: rp(:,:), rm(:,:), rr(:,:), r_norm(:,:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp), allocatable :: up(:,:), um(:,:), eigp(:,:), eigm(:,:), bp(:,:), bm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: epmat(:,:), emmat(:,:), smat(:,:)
!
!   auxiliary quantities for the helmich-paris algorithm
!
    real(dp), allocatable :: u_1(:,:), sv_1(:), vt_1(:,:), u_2(:,:), sv_2(:), vt_2(:,:), &
                             ept(:,:), emt(:,:), cmat(:,:), xpt(:,:), xmt(:,:), scr(:,:)
    real(dp), allocatable :: work_svd(:)
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
    lda2    = 2 * lda
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
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda2,lda2), a_copy(lda2,lda2), s_red(lda2,lda2), s_copy(lda2,lda2), &
      e_red(lda2), epmat(lda,lda), emmat(lda,lda), smat(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigp(n,n_max), eigm(n,n_max), bp(n,n_max), bm(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate additional memory for the helmich-paris algorithm:
!
    if (i_alg.eq.1) then
      allocate (u_1(lda,lda), sv_1(lda), vt_1(lda,lda), u_2(lda,lda), sv_2(lda), vt_2(lda,lda), &
                ept(lda,lda), emt(lda,lda), cmat(lda,lda), xpt(lda,lda), xmt(lda,lda), scr(lda,lda), stat = istat)
      call check_mem(istat)
!
      call dgesvd('a','a',lda,lda,u_1,lda,sv_1,u_1,lda,vt_1,lda,lw_svd,-1,info)
      lwork_svd = int(lw_svd(1))
      allocate (work_svd(lwork_svd))
    end if
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag  = zero
    t_ortho = zero
    t_mv    = zero
    t_tot   = zero
    vp      = zero
    vm      = zero
    bvp     = zero
    bvm     = zero
    lvp     = zero
    lvm     = zero
    a_red   = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max
      vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
      vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
    end do
!
    call ortho_cd(n,n_max,vp,growth,ok)
    call ortho_cd(n,n_max,vm,growth,ok)
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
      call apbmul(n,n_act,vp(1,i_beg),lvp(1,i_beg))
      call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
      call spdmul(n,n_act,vp(1,i_beg),bvm(1,i_beg))
      call smdmul(n,n_act,vm(1,i_beg),bvp(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,ldu,n,one,vp,n,lvp,n,zero,epmat,lda)
      call dgemm('t','n',ldu,ldu,n,one,vm,n,lvm,n,zero,emmat,lda)
      call dgemm('t','n',ldu,ldu,n,one,vm,n,bvm,n,zero,smat,lda)
!
      a_red = zero
      s_red = zero
      a_red(1:ldu,1:ldu)             = epmat(1:ldu,1:ldu)
      a_red(ldu+1:2*ldu,ldu+1:2*ldu) = emmat(1:ldu,1:ldu)
      s_red(1:ldu,ldu+1:2*ldu)       = transpose(smat(1:ldu,1:ldu))
      s_red(ldu+1:2*ldu,1:ldu)       = smat(1:ldu,1:ldu)
!
      call get_time(t1)
      if (i_alg.eq.0) then
!
!       default algorithm: solve the 2n-dimensional inverse problem.
!
        call dsygv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,work,lwork,info)
!
!       extract the eigenvalues and compute the ritz approximation to the
!       eigenvectors 
!
        do i_eig = 1, n_max
          eig(i_eig)      = one/e_red(2*ldu - i_eig + 1)
          up(1:ldu,i_eig) = s_red(1:ldu,2*ldu - i_eig + 1)
          um(1:ldu,i_eig) = s_red(ldu+1:2*ldu,2*ldu - i_eig + 1)
        end do
!
      else if (i_alg.eq.1) then
!     
!       helmich-paris algorithm
!
        scr = smat
!
!       compute the svd of smat
!
        call dgesvd('a','a',ldu,ldu,scr,lda,sv_1,u_1,lda,vt_1,lda,work_svd,lwork_svd,info)
!
!       divide the columns of u_1 and the rows of v_1 by the square root of the singular values
!
        do i = 1, ldu
          u_1(:,i) = u_1(:,i) / sqrt(sv_1(i))
          vt_1(i,:) = vt_1(i,:) / sqrt(sv_1(i))
        end do
!
!       for the projected ep and em matrices:
!
        call dgemm('n','t',ldu,ldu,ldu,one,epmat,lda,vt_1,lda,zero,scr,lda)
        call dgemm('n','n',ldu,ldu,ldu,one,vt_1,lda,scr,lda,zero,ept,lda)
        call dgemm('n','n',ldu,ldu,ldu,one,emmat,lda,u_1,lda,zero,scr,lda)
        call dgemm('t','n',ldu,ldu,ldu,one,u_1,lda,scr,lda,zero,emt,lda)
!
!       compute their cholesky decomposition:
!
        call dpotrf('l',ldu,ept,lda,info)
        call dpotrf('l',ldu,emt,lda,info)
!
!       assemble c = (l^-)^t l^+
!
        do j = 1, ldu
          do i = 1, ldu 
            cmat(i,j) = zero
            do k = max(i,j), ldu
              cmat(i,j) = cmat(i,j) + emt(k,i) * ept(k,j)
            end do
          end do
        end do
!
!       compute its svd:
!
        call dgesvd('a','a',ldu,ldu,cmat,lda,sv_2,u_2,lda,vt_2,lda,work_svd,lwork_svd,info)
!
!       assemble xpt and xmt
!
        do i = 2, ldu
          do j = 1, i-1
            ept(j,i) = zero
            emt(j,i) = zero
          end do
        end do
        call dgemm('n','n',ldu,ldu,ldu,one,emt,lda,u_2,lda,zero,scr,lda)
        call dgemm('t','n',ldu,ldu,ldu,one,vt_1,lda,scr,lda,zero,xpt,lda)
        call dgemm('n','t',ldu,ldu,ldu,one,ept,lda,vt_2,lda,zero,scr,lda)
        call dgemm('n','n',ldu,ldu,ldu,one,u_1,lda,scr,lda,zero,xmt,lda)
!
!       normalize and gather the eigenvalues and eigenvectors:
!
        do i_eig = 1, n_max
          eig(i_eig)  = sv_2(ldu - i_eig + 1)
          up(1:ldu,i_eig) = xpt(1:ldu,ldu - i_eig + 1)/(sqrt(two) * sv_2(ldu - i_eig + 1))
          um(1:ldu,i_eig) = xmt(1:ldu,ldu - i_eig + 1)/(sqrt(two) * sv_2(ldu - i_eig + 1))
        end do
!
!       gather the eigenvalues
!
      end if
      call get_time(t2)
      t_diag = t_diag + t2 - t1
      call dgemm('n','n',n,n_max,ldu,one,vp,n,up,lda,zero,eigp,n)
      call dgemm('n','n',n,n_max,ldu,one,vm,n,um,lda,zero,eigm,n)
!
      do i_eig = 1, n_max
        evec(1:n,i_eig)    = eigp(:,i_eig) + eigm(:,i_eig)
        evec(n+1:n2,i_eig) = eigp(:,i_eig) - eigm(:,i_eig)
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,lvp,n,up,lda,zero,rp,n)
      call dgemm('n','n',n,n_max,ldu,one,lvm,n,um,lda,zero,rm,n)
      call dgemm('n','n',n,n_max,ldu,one,bvp,n,um,lda,zero,bp,n)
      call dgemm('n','n',n,n_max,ldu,one,bvm,n,up,lda,zero,bm,n)
      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),bp(:,i_eig),1,rp(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bm(:,i_eig),1,rm(:,i_eig),1)
        r_norm(1,i_eig) = dnrm2(n,rp(:,i_eig),1)/sqrtn + dnrm2(n,rm(:,i_eig),1)/sqrtn
        r_norm(2,i_eig) = maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig)))
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
          write(6,1040) it, i_eig, eig(i_eig), r_norm(:,i_eig), done(i_eig)
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
        call lrprec(n,n_act,eig(ind),rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call ortho_vs_x(n,ldu,n_act,vp,vp(1,i_beg),xx,xx)
        call ortho_vs_x(n,ldu,n_act,vm,vm(1,i_beg),xx,xx)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!
        n_act = n_max 
        vp = zero
        vm = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
          vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
        end do
!
        call ortho_cd(n,n_max,vp,growth,ok)
        call ortho_cd(n,n_max,vm,growth,ok)
!
        lvp   = zero
        lvm   = zero
        bvp   = zero
        bvm   = zero
        a_red = zero
        s_red = zero
        epmat = zero
        emmat = zero
        smat  = zero
!
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
    1000 format(t3,'timings for caslr (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
    deallocate(work,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm,a_red,a_copy,s_red,s_copy,e_red, &
               epmat,emmat,smat,up,um,eigp,eigm,bp,bm)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_driver
!
  subroutine caslr_eff_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                              apbmul,ambmul,spdmul,smdmul,lrprec,eig,evec,ok)
!
!   main driver for the efficient solution to the following generalized 
!   eigenvalue problem:
!
!   / A  B \ / Y \     /  S  D \ / Y \ 
!   |      | |   | = w |       | |   |
!   \ B  A / \ Z /     \ -D -S / \ Z /
!
!   where A, B, S are symmetric matrices and D is antysimmetric.
!
!   If (w, Y, Z) are a solution, then (-w, Z, Y) is also a solution. 
!   Following J. Chem. Phys., 118, 522 (2003), we enforce this property in 
!   the iterative procedure by expanding the eigenvector as
!
!   (Y, Z) = (b+, b+) + (b-, - b-)
!
!   This routine solves the associate problem 
!
!   /  S  D \ / Y \   1 / A  B \ / Y \ 
!   |       | |   | = - |      | |   |
!   \ -D -S / \ Z /   w \ B  A / \ Z /
!
!   using the casida matrix, which is symmetric and positive definite, as 
!   the metric. This allows us to use expansion vectors that are orthogonal
!   with respect to the dot product defined by the metric, which in turn
!   results in a Rayleigh-Ritz procedure that requires the solution of a 
!   symmetric standard eigenvalue problem
!
!   / 0  s^T \ / u+ \   1 / u+ \
!   |        | |    | = - |    |
!   \ S   0  / \ u- /   w \ u- /
!
!   which can be reduced to a half-sized eigenvalue problem
!
!   s^T s u+ = (1/w)^2 u+
!   u- = 1/w Su+
!
!   input variables
!   ===============
!
!   verbose:  logical, whether to print various information at each 
!             iteration (eigenvalues, residuals...).
!
!   n:        integer, size of the A, B, S, D matrices
!
!   n2:       integer, n2 = 2*n, size of the generalized eigenvalue
!             problem
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
!   maxdav:   maximum size of the expansion subspace
!
!   apbmul:   external subroutine to apply (A+B) to a vector
!
!   ambmul:   external subroutine to apply (A-B) to a vector
!
!   spdmul:   external subroutine to apply (S+D) to a vector
!
!   smdmul:   external subroutine to apply (S-D) to a vector
!
!   lrprec:   external subroutine that applies a preconditioner.
!
!   output variables:
!   =================
!
!   eig:      double precision real array of size n_max. 
!             if ok is true, contains the positive eigenvalues
!             of the generalized problem in ascending order.
!
!   evec:     double precision real array of size (n2,n_max).
!             in input, a guess for the eigenvectors.
!             if ok is true, in output the computed eigenvectors.
!
!   ok:       logical, true if caslr_eff_driver converged.
!
    use utils
    implicit none
    logical, intent(in)                          :: verbose
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n_max),    intent(inout) :: eig
    real(dp), dimension(n2,n_max), intent(inout) :: evec
    logical,                       intent(inout) :: ok
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
!
!   actual expansion space size and total dimension
!
    integer               :: dim_dav, lda, lda2
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
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    real(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    real(dp), allocatable :: rp(:,:), rm(:,:), rr(:,:), r_norm(:,:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp), allocatable :: up(:,:), um(:,:), eigp(:,:), eigm(:,:), bp(:,:), bm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: smat(:,:)
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
    lda2    = 2 * lda
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
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (s_red(lda,lda), s_copy(lda,lda), e_red(lda2), smat(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigp(n,n_max), eigm(n,n_max), bp(n,n_max), bm(n,n_max), stat = istat)
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
    t_diag  = zero
    t_ortho = zero
    t_mv    = zero
    t_tot   = zero
    vp      = zero
    vm      = zero
    bvp     = zero
    bvm     = zero
    lvp     = zero
    lvm     = zero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max
      vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
      vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
    end do
!
!   orthogonalize the expansion space to the metric.
!
    call apbmul(n,n_max,vp,lvp)
    call b_ortho(n,n_max,vp,lvp)
    call ambmul(n,n_max,vm,lvm)
    call b_ortho(n,n_max,vm,lvm)
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
!     perform this iteration's matrix-vector multiplications:
!
      call get_time(t1)
      call spdmul(n,n_act,vp(1,i_beg),bvm(1,i_beg))
      call smdmul(n,n_act,vm(1,i_beg),bvp(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,ldu,n,one,vm,n,bvm,n,zero,smat,lda)
!
!     save s, and assemble s^t s:
!
      s_red  = smat
      s_copy = zero
      call dgemm('t','n',ldu,ldu,ldu,one,s_red,lda,s_red,lda,zero,s_copy,lda)
!
!     diagonalize s^t s
!
      call get_time(t1)
      call dsyev('v','u',ldu,s_copy,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      do i_eig = 1, n_max
        eig(i_eig)      = sqrt(e_red(ldu - i_eig + 1))
        up(1:ldu,i_eig) = s_copy(1:ldu,ldu - i_eig + 1)
      end do
!
!     compute the u_- eigenvectors:
!
      call dgemm('n','n',ldu,n_max,ldu,one,s_red,lda,up,lda,zero,um,lda)
      do i_eig = 1, n_max
        um(1:ldu,i_eig) = um(1:ldu,i_eig)/eig(i_eig)
      end do
!
!     asemble the symmetric and antysimmetric combinations (Y+Z) and (Y-Z)
!
      call dgemm('n','n',n,n_max,ldu,one,vp,n,up,lda,zero,eigp,n)
      call dgemm('n','n',n,n_max,ldu,one,vm,n,um,lda,zero,eigm,n)
!
!     assemble the current approximation to the eigenvectors
!
      do i_eig = 1, n_max
        evec(1:n,i_eig)    = eigp(:,i_eig) + eigm(:,i_eig)
        evec(n+1:n2,i_eig) = eigp(:,i_eig) - eigm(:,i_eig)
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ldu,one,bvp,n,um,lda,zero,rp,n)
      call dgemm('n','n',n,n_max,ldu,one,bvm,n,up,lda,zero,rm,n)
      call dgemm('n','n',n,n_max,ldu,one,lvp,n,up,lda,zero,bp,n)
      call dgemm('n','n',n,n_max,ldu,one,lvm,n,um,lda,zero,bm,n)
      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),bp(:,i_eig),1,rp(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bm(:,i_eig),1,rm(:,i_eig),1)
        r_norm(1,i_eig) = (dnrm2(n,rp(:,i_eig),1) + dnrm2(n,rm(:,i_eig),1))/(eig(i_eig)*sqrt(two)*sqrtn)
        r_norm(2,i_eig) = (maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig))))/(sqrt(two)*eig(i_eig))
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
          write(6,1040) it, i_eig, one/eig(i_eig), r_norm(:,i_eig), done(i_eig)
        end do
        write(6,*) 
      end if
!
      if (all(done(1:n_targ))) then
        ok = .true.
        do i_eig = 1, n_targ
          eig(i_eig) = one/eig(i_eig)
        end do
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
        call lrprec(n,n_act,eig(ind),rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call b_ortho_vs_x(n,ldu,n_act,vp,lvp,vp(1,i_beg))
        call apbmul(n,n_act,vp(1,i_beg),lvp(1,i_beg))
        call b_ortho(n,n_act,vp(1,i_beg),lvp(1,i_beg))
        call b_ortho_vs_x(n,ldu,n_act,vm,lvm,vm(1,i_beg))
        call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call b_ortho(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
!
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       initialize indexes back to their starting values 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!        
        n_act = n_max 
        vp = zero
        vm = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
          vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
        end do
!
        lvp   = zero
        lvm   = zero
!
        call apbmul(n,n_max,vp,lvp)
        call b_ortho(n,n_max,vp,lvp)
        call ambmul(n,n_max,vm,lvm)
        call b_ortho(n,n_max,vm,lvm)
!
        bvp   = zero
        bvm   = zero
        s_red = zero
        smat  = zero
!        
      end if
!      
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!      
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
    1000 format(t3,'timings for caslr_eff (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
    deallocate(work,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm,s_red,s_copy,e_red,smat,up,um,eigp,eigm,bp,bm)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_eff_driver
!
  subroutine davidson_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav,&
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
        call ortho_vs_x(n,ldu,n_act,space,space(1,i_beg),xx,xx)
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
!   apply the b matrix and b-orthogonalize the guess:
!
    call bvec(n,n_max,space,bspace)
    call b_ortho(n,n_max,space,bspace)
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
!     call bvec(n,n_act,space(1,i_beg+n_rst),bspace(1,i_beg+n_rst))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrices 
!
      call dgemm('t','n',ldu,n_act,n,one,space,n,aspace(1,i_beg+n_rst),n,zero,a_red(1,i_beg+n_rst),lda)
!     call dgemm('t','n',ldu,n_act,n,one,space,n,bspace(1,i_beg+n_rst),n,zero,s_red(1,i_beg+n_rst),lda)
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
!     s_copy = s_red
!
!     diagonalize the reduced matrix
!
      call get_time(t1)
!     call dsygv(1,'v','u',ldu,a_copy,lda,s_copy,lda,e_red,work,lwork,info)
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
        call b_ortho_vs_x(n,ldu,n_act,space,bspace,space(1,i_beg))
        call bvec(n,n_act,space(1,i_beg),bspace(1,i_beg))
        call b_ortho(n,n_act,space(1,i_beg),bspace(1,i_beg))
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
        call b_ortho(n,n_max,space,bspace)
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
  subroutine nonsym_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav,shift,&
                            matvec,matvec_l,precnd,eig,evec_r,evec_l,both,ok)
    logical,                        intent(in)    :: verbose
    integer,                        intent(in)    :: n, n_targ, n_max
    integer,                        intent(in)    :: max_iter, max_dav
    real(dp),                       intent(in)    :: tol, shift
    real(dp), dimension(n_max),     intent(inout) :: eig
    real(dp), dimension(n,n_max),   intent(inout) :: evec_r, evec_l
    logical,                        intent(inout) :: both, ok
    external                                      :: matvec, matvec_l, precnd
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
    real(dp)              :: sqrtn, tol_rms, tol_max, tol_im
    real(dp)              :: xx(1)
!
!   arrays to control convergence and orthoginalization
!
    logical, allocatable  :: done(:)
!
!   expansion spaces, residuals and their norms
!
    real(dp), allocatable :: space_r(:,:), space_l(:,:), aspace_r(:,:), aspace_l(:,:), &
                              r_r(:,:), r_l(:,:), r_norm_r(:,:), r_norm_l(:,:)
!
!   subspace matrix, eigenvalues and real and imaginary parts of the eigenvalues
!
    real(dp), allocatable :: a_red_r(:,:), a_red_l(:,:), e_red_re(:), e_red_im(:), &
                             evec_red_r(:,:), evec_red_l(:,:), copy_r(:,:), copy_l(:,:), copy_eig(:), overlap(:,:)
!
!   restarting variables
!
    logical               :: restart
    integer               :: n_rst
!
!   variables for orthogonalization
!
    real(dp),    allocatable  :: xu(:,:), qr(:,:), v(:,:), metric(:,:), msave(:,:), vy(:,:)
    integer(dp), allocatable  :: ipiv(:)
    logical                   :: done_lu(2)
    real(dp)                  :: vy_norm, tol_ortho_lu, alpha, unorm, shift_lu
    logical                   :: not_orthogonal, use_qr, micro_done, ok_lu
    integer                   :: k, i, j, max_orth, max_GS, qr_dim, it_micro, maxit_lu
!
    integer                   :: verbosity
! 
!   variables for the sort
!
    integer                   :: max_idx(1), fin, iter
    real(dp)                  :: diff(n_max), temp, t_sort(2)
    logical                   :: found_im, found_er, no_match
    logical, allocatable      :: mask_sort(:)
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2
    external              :: dnrm2, dgemm, dcopy, dgeev, dsyev
!
!   computing actual size of the expansion space, checking that
!   the input makes sense.
!
    dim_dav =  max(min_dav,max_dav)
    lda = dim_dav*n_max 
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max) 
    allocate (work(lwork), tau(lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for for expansion space, the corresponding
!   matrix-multiplied vectors and the residuals
!
    allocate (space_r(n,lda), space_l(n,lda), aspace_r(n,lda), aspace_l(n,lda), r_r(n,n_max), &
      r_l(n,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for convergency check
!
    allocate (done(n_max), r_norm_r(2,n_max), r_norm_l(2,n_max), stat=istat)
!
!   allocate memory for the reduced matrix, its eigenvalues with real &
!   imaginary parts, and its left & right eigenvectors 
!
    allocate (a_red_r(lda,lda), a_red_l(lda,lda), e_red_re(2*lda), e_red_im(2*lda), evec_red_r(lda,lda), &
              evec_red_l(lda,lda), copy_r(lda,lda), copy_l(lda,lda), copy_eig(2*lda), overlap(n_max,n_max), stat =istat)
    call check_mem(istat)
!
!   allocate space for orthogonalization routines
!
    allocate (xu(lda,n_max))
    allocate (vy(lda,n_max))
    allocate (qr(n,lda))
!   
!   allocate mask array for sorting routine
!
    allocate (mask_sort(lda))
!
!   set the tolerance and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
    tol_ortho_lu = 1.d-14
    tol_im  = 1.d-12
!
!   clean out various quantities
!
    t_diag    = zero
    t_ortho   = zero
    t_mv      = zero
    t_tot     = zero
    t_sort    = zero
    space_r   = zero
    space_l   = zero
    aspace_r  = zero
    aspace_l  = zero
    a_red_r   = zero
    a_red_l   = zero
    e_red_re  = zero
    e_red_im  = zero
    copy_r    = zero
    copy_l    = zero
    r_norm_r  = zero
    r_norm_l  = zero
    done      = .false.
    done_lu   = .true. 
    use_qr    = .false.
    verbosity = 0
    ok        = .false.
!
    call get_time(t_tot)
!
!   check weather we have a guess for the eigenvectors in evec, and
!   weather it is orthonormal.
!   if evec is zero, create a random guess
!
    call check_guess(n,n_max,evec_r)
    call check_guess(n,n_max,evec_l)
!
!   move guess into the expansion spaces
!
    call dcopy(n*n_max,evec_r,1,space_r,1)
    call dcopy(n*n_max,evec_l,1,space_l,1)
!
!   initialize the number of active vectors and the associated indices.
!
    n_act = n_max
    i_beg = 1
    ind   = 1
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
!   main loop
!
    1030 format(t5,'Davidson-nonsym iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
! 
    if (verbose) write(6,1030) tol
!   
    n_rst = 0
    do it = 1, max_iter
!
!     header
!
      if (verbosity.ge.1) then
        print *
        1000 format(20("="),/,t3,'Iteration',i3,/,20("="))
        print 1000, it
        print*
      end if
!
!     update the size of the expansion space.
!
      ldu   = ldu + n_act
!
!     perform this iteration's matrix-vector multiplications for both 
!     right and left spaces
!
      call get_time(t1)
      call matvec(n,n_act,space_r(1,i_beg),aspace_r(1,i_beg))
      if (both) call matvec_l(n,n_act,space_l(1,i_beg),aspace_l(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 -t1
!
      if (verbosity.ge.2) then
        1200 format(t3, 'ldu:     ',i4,/, &
                    t3, 'i_beg:   ', i4,/,&
                    t3, 'm_dim:   ', i4,/,&
                    t3, 'n_rst:   ', i4,/,&
                    t3, 'restart: ', l4,/,&
                    t3, 'n_act:   ', i4,/)
        write(6,1200) ldu, i_beg, m_dim, n_rst, restart, n_act
        print *, "space r+l"
        call printMatrix(n,ldu,space_r,n)
        print *
        call printMatrix(n,ldu,space_l,n)
        print *
        print *, "aspace r+l"
        call printMatrix(n,ldu,aspace_r,n)
        print *
        call printMatrix(n,ldu,aspace_l,n)
        print *
      end if
!
!     get the reduced matrix
!
      a_red_r = zero
      if (both) then
        call dgemm('t','n',ldu,ldu,n,one,space_l,n,aspace_r,n,zero,a_red_r,lda)
      else
        call dgemm('t','n',ldu,ldu,n,one,space_r,n,aspace_r,n,zero,a_red_r,lda)
      end if
!
      if (verbosity.ge.1) then
        print *
        print *, "print reduced space a:"
        call printMatrix(ldu,ldu,a_red_r,lda)
        print *
      end if
!
!     diagonalize the reduced matrix
!
      call get_time(t1)
      !if (both) then
      call dgeev('v','v',ldu,a_red_r,lda,e_red_re,e_red_im,evec_red_l,lda,evec_red_r,lda,work,lwork,info)
       !else
      !  evec_red_r = a_red_r
      !  call dsyev('v','l',ldu,evec_red_r,lda,e_red_re,work,lwork,info)
      !end if
!
!     debug symmetric case
!
      call get_time(t2)
!
      t_diag = t_diag + t2 - t1
! 
!     check if diagonalization terminated with info = 0 and check if 
!     the eigenvalues have a complex contribution
!
      if (info.ne.0) then
        print *, "diagonalization of reduced space failed."
        stop
      end if
!
!     sort eigenvalues and eigenvectors in decreasing order in range n_targ 
!
      if (verbosity.ge.2) then
        print *, "before sorting"
        print *, "evec r + l"
        call printMatrix(ldu,ldu,evec_red_r,lda)
        print *
        call printMatrix(ldu,ldu,evec_red_l,lda)
        print *
        print*, "eig"
        call printVector(e_red_re,ldu)
        print *
      end if
!
!     test sort function with test set
!
      call get_time(t1)
      call sort_eigenpairs(ldu,ldu,e_red_re,e_red_im,evec_red_r,evec_red_l,n_max,lda,.true.,tol_im)
      call get_time(t2)
      t_sort = t_sort + t2 - t1
!
!     double check for complex contributions in the n_max sought eigenvalues
!
      found_im = .false.
      do j = 1, n_max
        if (e_red_im(j).gt.tol_im) found_im = .true.
      end do
!
      if (found_im.and.verbose) then
        print *
        print *, "complex contribution in sought eigenvalues"
        print*
      end if
!
      if (verbosity.ge.2) then
        print *
        print *, "eigenvalues after sort real and imaginary"
        call printVector(e_red_re, ldu)
        print *
        call printVector(e_red_im, ldu)
        print *
      end if
!
!
      if (verbosity.ge.2) then
        print *, "evec r + l after sort + eigenvecs"
        call printMatrix(ldu,ldu,evec_red_r,lda)
        print *
        call printMatrix(ldu,ldu,evec_red_l,lda)
        print *
        call printVector(e_red_re,ldu)
        print *
      end if
!
!     compute overlap of old and new eigenvectors in the dimension of the old eigenvectors
!     to ensure correct sorting by checking if largest absolute value of column is on the
!     diagonal
!  
      if (it.ne.1 .and. .not. restart) then  
        mask_sort = .true.
        no_match  = .true.
        iter      = 0
        do while (no_match)
          iter = iter +1
          if (iter .gt. 10*n_max) then
            print *, "too many iterations in shifting away eigenvalues we dont like."
            stop
          end if
!
          call dgemm('t','n',n_max,n_max,ldu,one,copy_r,lda,evec_red_r,lda,zero,overlap,n_max)
!
          found_er = .false.
          do j = 1, n_max
            max_idx = maxloc(abs(overlap(:,j)))
            if (max_idx(1).ne.j) then
              found_er = .true.
              !call printMatrix(n_max,n_max,overlap,n_max)
              !stop
            end if
          end do
!
          call dgemm('t','n',n_max,n_max,ldu,one,copy_l,lda,evec_red_l,lda,zero,overlap,n_max)
!
          do j = 1, n_max
            max_idx = maxloc(abs(overlap(:,j)))
            if (max_idx(1).ne.j) then
              found_er = .true.
              !call printMatrix(n_max,n_max,overlap,n_max)
              !stop
            end if
          end do
!
          if (found_er) then
!
!           if error encounterd in the overlap of old and new eigenvectors, move errorneous 
!           eigenpair to the end of the array, mask it and sort again
!
!           identify difference of every old and new eigenvalue in range n_max and store 
!           lowest difference to not shift-away sought eigenvalues, which were just returned
!           in different order after the diagonalization, but are correct ones.
!
            diff    = 0
!
            do j = 1, n_max
              do k = 1, ldu
                temp = abs(e_red_re(j) - copy_eig(k)) 
                if (temp.lt.diff(j) .or. k.eq.1) diff(j) = temp
              end do
            end do
!
!           get index of highest difference and shift eigenpair to last available entry of array 
!           in range ldu and mask it
!
            max_idx = maxloc(diff)
            fin = ldu
            do j =1, ldu
              if (.not. mask_sort(j) .and. fin.gt.1) then
                fin = fin - 1
              else if (fin.le.1) then
                print *, "eigenvectors dont match although handling was tried."
                stop
              end if 
            end do
!
            call get_time(t1)
            call swap_eigenpairs(max_idx(1),fin,ldu,ldu,e_red_re,e_red_im,evec_red_r,evec_red_l,lda)
            call get_time(t2)
            t_sort = t_sort + t2 - t1
            mask_sort(fin) = .false.
!
            call get_time(t1)
            call sort_eigenpairs(ldu,ldu,e_red_re,e_red_im,evec_red_r,evec_red_l,n_max,lda,.true.,tol_im,mask_sort)
            call get_time(t2)
            t_sort = t_sort + t2 - t1
!
!           
            if (verbosity.gt.2) then
              print*
              print *, "---- Information ----"
              print *, "handled inconsistance in old and current eigenvectors"
              print *
            end if
          else
            no_match = .false.      
          end if
        end do
      end if
!
!     copy and save the new eigenvectors for the next iteration
!
      restart = .false.
      copy_r   = evec_red_r
      copy_l   = evec_red_l
      copy_eig = e_red_re
!
!     extract the eigenvalues and compute the ritz approximation to the 
!     eigenvectors
!
      eig = e_red_re(1:n_max)
!
      call dgemm('n','n',n,n_max,ldu,one,space_r,n,evec_red_r,lda,zero,evec_r,n)
      if (both) call dgemm('n','n',n,n_max,ldu,one,space_l,n,evec_red_l,lda,zero,evec_l,n)
!
      if (verbosity.ge.2) then
        print *
        print *, "extracted eigenvals"
        print *
        call printVector(eig, n_max)
        print *
!
        print *, "ritz"
        call printMatrix(n,n_max,evec_r,n)
        print *
        call printMatrix(n,n_max,evec_l,n)
        print *
      end if
!
!     compute the residuals, and their rms and sup norms
!
      call dgemm('n','n',n,n_max,ldu,one,aspace_r,n,evec_red_r,lda,zero,r_r,n) 
      if (both) call dgemm('n','n',n,n_max,ldu,one,aspace_l,n,evec_red_l,lda,zero,r_l,n) 
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n,-eig(i_eig),evec_r(:,i_eig),1,r_r(:,i_eig),1)
        r_norm_r(1,i_eig) = dnrm2(n,r_r(:,i_eig),1)/sqrtn
        r_norm_r(2,i_eig) = maxval(abs(r_r(:,i_eig)))
        if (both) then
          call daxpy(n,-eig(i_eig),evec_l(:,i_eig),1,r_l(:,i_eig),1)
          r_norm_l(1,i_eig) = dnrm2(n,r_l(:,i_eig),1)/sqrtn
          r_norm_l(2,i_eig) = maxval(abs(r_l(:,i_eig)))
        end if 
!
      end do
!
      if (verbosity.ge.2) then
        print*
        print*, "residual right and left"
        call printMatrix(n,n_max,r_r,n)
        print*
        call printMatrix(n,n_max,r_l,n)
        print*
      end if
!
!     check convergency. look the first contiguous converged eigenvalues
!     by setting the logical array "done" to true
!
      do i_eig = 1, n_targ
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm_r(1,i_eig).lt.tol_rms .and.&
                          r_norm_r(2,i_eig).lt.tol_max .and.&
                          r_norm_l(1,i_eig).lt.tol_rms .and.& 
                          r_norm_l(2,i_eig).lt.tol_max .and.&    
                          it.gt.1                               
        if (.not.done(i_eig)) then
          done(i_eig+1:n_max) = .false.
          exit
        end if
      end do
!
!     print some information
!
      if (verbose) then
        if (verbosity.ge.1) print *, "interim info"
        do i_eig = 1, n_targ
          write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm_r(:,i_eig), done(i_eig)
        end do
        write(6,*)
      end if 
!
      if (all(done(1:n_targ))) then
        ok = .true.
        if (verbose) print *, "converged =)"
        exit
      end if
!     
!     check weather an update is required.
!     if not, perform a davidson restart
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
        call precnd(n,n_act,-eig(ind),r_r(1,ind),space_r(1,i_beg))
        if (both) call precnd(n,n_act,-eig(ind),r_l(1,ind),space_l(1,i_beg))
!
        if (verbosity.ge.2) then
          print*
          print *, "precondition:"
          print *
          call printMatrix(n,n_act,space_r(1,i_beg),n)
          print*
          print*
          call printMatrix(n,n_act,space_l(1,i_beg),n)
          print*
        end if
!
!       orthogonalize the new vectors to the existing ones of the respective other
!       space and orthogonalize set of new vectors among each other
!
        call get_time(t1)
!
!       Gram-Schmit orthogonalization of residual to the respective subspace 
!
        if (both) then
          call biortho_vs_x(n,ldu,n_act,space_l,space_r,space_l(1,i_beg),space_r(1,i_beg))
        else
          call ortho_vs_x(n,ldu,n_act,space_r,space_r(1,i_beg),xx,xx)
        end if
!
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
        !do i = 0, n_max-1
        !  space_l(:,i_beg+i) = space_l(:,i_beg+i)/dnrm2(n,space_l(:,i_beg+i),1)
        !  space_r(:,i_beg+i) = space_r(:,i_beg+i)/dnrm2(n,space_r(:,i_beg+i),1)
        !end do
!
        if (verbosity.ge.2) then
          call printMatrix(n,ldu+n_max,space_r,n)
          print*
!
          print * 
          call printMatrix(n,ldu+n_max,space_l,n)
          print*
          print *
        end if
!
!
!       normalize columns 
!
        if (verbosity.gt.2) then
          print *, "normalized space r"
          call printMatrix(n,ldu+n_max,space_r,n)
          print *, "normalized space l"
          call printMatrix(n,ldu+n_max,space_l,n)
        end if
!
      else 
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        n_act   = n_max
        space_r = zero
        space_l = zero
!
!       put current eigenvectors into the first position of tne
!       expansion space
!
        call dcopy(n_max*n,evec_r,1,space_r,1)
        call dcopy(n_max*n,evec_l,1,space_l,1)
!
        call get_time(t1)
        if (both) call svd_biortho(n,n_act,space_r,space_l)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
        aspace_r = zero
        aspace_l = zero
        a_red_r  = zero
        a_red_l  = zero
        e_red_re  = zero
        e_red_im  = zero
!
!       initialize indexes back to their starting values
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
!
!       counting how many matvec we can skip at the next
!       iteration
!
        do i_eig = 1, n_targ
          if (done(i_eig)) then
            n_frozen = n_frozen +1
          else
            exit
          end if
        end do
        restart = .true.
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
    end do
!
!   if required, print timings
!
    1100 format(t3,'timings for davidson (cpu/wall): ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'  sorting:                       ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1100) t_mv, t_diag, t_ortho, t_sort, t_tot
!      
!   deallocate memory
!
    deallocate(work, tau, space_r, space_l, aspace_r, aspace_l, r_r, r_l, done, r_norm_r, r_norm_l)
    deallocate(a_red_r, a_red_l, e_red_re, e_red_im, evec_red_r, evec_red_l, copy_r, copy_l)
    deallocate(xu, vy, qr)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
!
  end subroutine nonsym_driver
!
  subroutine sort_eigenpairs(n,m,wr,wl,vr,vl,n_want,ldv,ignore,thresh,mask_in)
!
!   sort m real & imaginary eigenvalues and right & left eigenvectors of length n 
!   in decreasing order according to the real eigenvalues in the rang.e of n_want
! 
    implicit none
    integer,  intent(in)      :: n, m, ldv, n_want
    real(dp), intent(inout)   :: wr(m), wl(m), vr(ldv,m), vl(ldv,m)
    real(dp), intent(in)      :: thresh
    logical,  intent(in)      :: ignore
    logical,  optional        :: mask_in(m)
!   
!   local variables
!
    real(dp)                  :: w, v(ldv)
    integer                   :: i, j, idx, min_idx(1), fin
    logical                   :: mask(m)
!
    real(dp)                  :: dnrm2
!
!   define initial mask
!
    if (present(mask_in)) then
      mask = mask_in
    else
      mask = .true.
    end if
!
    do i = 1, n_want
! 
!     identify minimal value and mask first position for next iteration
!
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
        call swap_eigenpairs(fin,idx,n,m,wr,wl,vr,vl,ldv)
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
      call swap_eigenpairs(i,idx,n,m,wr,wl,vr,vl,ldv)
! 
!
    end do
!
  end subroutine sort_eigenpairs
!
  subroutine swap_eigenpairs(i,j,n,m,wr,wl,vr,vl,ldv)
!
!   swaps m real & immaginary eigenvalues and eigenvectors of length l of the
!   indices i and j with each other 
!
    implicit none
    integer, intent(in)       :: n, m, ldv, i, j
    real(dp), intent(inout)   :: wr(m), wl(m), vr(ldv,m), vl(ldv,m)
!
    real(dp)                  :: w, v(ldv)
!
    w       = wr(i)
    wr(i)   = wr(j)
    wr(j)   = w
!
    w       = wl(i)
    wl(i)   = wl(j)
    wl(j)   = w
!
    v         = vr(:,i)
    vr(:,i)   = vr(:,j)
    vr(:,j)   = v
!   
    v         = vl(:,i)
    vl(:,i)   = vl(:,j)
    vl(:,j)   = v
!
    return
  end subroutine 
!
  subroutine checkInfo(info, occasion)
!    
! check if info is zero and print error message if not
!
    implicit none
    integer              :: info, zero
    character(len=*)     :: occasion

    zero = 0
    if (info .NE. zero) then
      print *
      print *, '--- WARNING ---'
      print *, occasion
      print *, 'Process terminated with info not equal to 0'
    end if
!
  end subroutine checkInfo
!
  subroutine printMatrix(n,m,A,lda) 
!   
! print formatted matrix
!
    implicit none
    integer , intent(in)  :: n, m, lda
    real(dp), intent(in)  :: A(lda,lda)
!
    integer :: i, j
!
    do i = 1, n
      do j = 1, m
        write(*,'(E15.5)', advance='no') A(i,j)
        if (j .lt. m) then
          write(*, '(A)', advance='no') ' '
        end if
      end do
      print *
    end do
!
  end subroutine printMatrix
!
  subroutine printVector(vec, m)
!    
! print formatted vector
!
    implicit none
    integer,  intent(in)  :: m
    real(dp), intent(in)  :: vec(m)
    integer               :: i
!
    do i = 1, m
      write(*,'(F13.5)') vec(i)
    end do
!
  end subroutine printVector
!
  subroutine printPythonMat(n,m,a,lda)
!   
! print formatted matrix
!
    integer , intent(in)  :: n, m, lda
    real(dp), intent(in)  :: a(lda,m)
!
    integer :: i,j 
    character(len=4000)  :: str
    character(len=30) :: temp
!
    str = "mat = np.array(["
    do i = 1, n
      str = trim(str) // ' ' // trim("[")
      do j = 1, m
        write(temp, '(F14.10)') a(i,j)
        str = trim(str) //  ' ' //trim(temp) // trim(',')
      end do
      str = trim(str) // ' ' //  trim("],")
    end do
    str = trim(str) // ' ' // trim("])")
    print *, str
!
  end subroutine printPythonMat
!
  subroutine checkOrth1mat(mat1,n,m,mat1_name,thresh,not_orthogonal,verbose)
! check if the vectors within the matrix are orthogonal
!
    implicit none
    integer,              intent(in)      ::  m,n
    real(dp),             intent(in)      ::  thresh, mat1(n,m)
    character(len=*),     intent(in)      ::  mat1_name
    logical,              intent(inout)   ::  not_orthogonal
    logical,              intent(in)      ::  verbose
    real(dp)                              ::  dot_prod
    integer :: i,j 

    do i = 1, m
      do j = 1, m
        if (j .NE. i) then
          dot_prod  = abs(dot_product( mat1(:,j), mat1(:,i) ))
          if ( dot_prod .GT. thresh ) then
            not_orthogonal = .true.
            if (verbose) then 
              print *
              print *, '--- WARNING ---'
              print *, 'vector ', i, ' of matrix ', mat1_name, ' is not orthogonal to vector ', j, ' of matrix', mat1_name
              print *, 'Result of dot product: ', dot_prod
              print *
            end if
          end if
        end if
      end do
    end do

  end subroutine checkOrth1mat
!
  subroutine checkOrth2mat(mat1, m1, mat1_name, mat2, m2, mat2_name, n, thresh, not_orthogonal, verbose)
!    
! check if vectors of two matrices are orthogonal
!
    implicit none
    integer,              intent(in)      :: m1, m2, n
    real(dp),             intent(in)      :: thresh, mat1(n,m1), mat2(n,m2)
    character(len=*),     intent(in)      :: mat1_name, mat2_name
    real(dp)                              :: dot_prod
    integer                               :: i,j 
    logical,              intent(inout)   :: not_orthogonal
    logical,              intent(in)      :: verbose
!
    do i = 1, m1
      do j = 1, m2
        dot_prod = abs(dot_product( mat2(:,j), mat1(:,i) ))
        if ( dot_prod .GT. thresh ) then
          not_orthogonal = .true.
          if (verbose) then
            print *
            print *, '--- WARNING ---'
            print *, 'vector ', i, ' of matrix ', mat1_name, ' is not orthogonal to vector ', j, ' of matrix', mat2_name
            print *, 'Result of dot product:', dot_prod
            print *
            !call printVector(mat1(:,i),n)
            !print *
            !call printVector(mat2(:,j),n)
          end if 
        end if
      end do
    end do
    !if (not_orthogonal) stop 'Not ortho'
!
  end subroutine checkOrth2mat
!
! orthogonalization routines:
! ===========================
!
  subroutine ortho(n,m,u,w)
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
    logical,  parameter   :: use_svd = .false.
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
  subroutine ortho_cd(n,m,u,growth,ok)
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
!   as cholesky decomposition is not the most stable way of orthogonalizing
!   a set of vectors, the orthogonalization is refined iteratively. 
!   a conservative estimate of the orthogonalization error is used to 
!   assess convergence. 
!
!   this routine returns a growth factor, which can be used in (b_)ortho_vs_x
!   to estimate the orthogonality error introduced by ortho_cd.
!
!   while it is very unlikely to do so, this routine can fail. 
!   a logical flag is then set to false, so that the calling program can 
!   call a more robust orthogonalization routine without aborting.
!
!   arguments:
!   ==========
!
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(inout) :: u
    real(dp),                  intent(inout) :: growth
    logical,                   intent(inout) :: ok
!
!   local variables
!   ===============
!
    integer               :: it, it_micro
    real(dp)              :: error, dnrm2, alpha, unorm, shift
    real(dp)              :: rcond, l_norm, linv_norm
    logical               :: macro_done, micro_done
    real(dp), parameter   :: tol_ortho_cd = two * epsilon(one)
    integer,  parameter   :: maxit = 10
!
!   local scratch
!   =============
!
    real(dp), allocatable :: metric(:,:), msave(:,:)
!
!   external functions:
!   ===================
!
    external              :: dgemm, dgeqrf, dtrsm, dnrm2, dtrcon
!
!   get memory for the metric.
!
    allocate (metric(m,m), msave(m,m))
    metric = zero
    macro_done = .false.
!
!   assemble the metric
!
    it = 0
    growth = one
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
      msave = metric
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
            stop
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
!     we assume that the error on the orthogonality is of order k(l)^2 * eps,
!     where eps is the machine precision. 
!     the condition number k(l) is estimated by computing 
!
!     k(l) ||l|| ||l^-1||,
!
!     where the norm used is the following (see norm_estimate):
!
!     || A || = || D + O || <= || D ||_inf + || O ||_2
!
!     compute l^-1, using msave to store the inverse cholesky factor
!
      msave = metric
      call dtrtri('l','n',m,msave,m,info)
!
!     compute the norm of l, l^-1 and the condition number:
!
      l_norm    = norm_est(m,metric)
      linv_norm = norm_est(m,msave)
      rcond     = l_norm * linv_norm
!
!     in each iteration of ortho_cd, we apply l^-t to u, which introduces
!     a numerical error of order ||l^-1||. 
!     this error is saved in growth and used in ortho_vs_x to check how much
!     ortho_cd spoiled the previously computed orthogonality to x. 
!
      growth    = growth * linv_norm
!
!     orthogonalize u by applying l^(-t)
!    
      call dtrmm('r','l','t','n',n,m,one,msave,m,u,n)
!
!     check the error:
!
      error      = epsilon(one) * rcond*rcond
      macro_done = error .lt. tol_ortho_cd
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
  subroutine ortho_lu(n,m,u_l,u_r,ok)
  implicit none

! 
!   orthogonalize m vectors of length n using the Cholesky factorization
!   of their overlap.
!   this is done by metric = U_l^t U_r and then by computing its lu 
!   decomposition metric = L U. The orthogonal vectors are obtained then
!   by solving the triangular linear systems
! 
!     U_l(ortho) L^t = U_l
!     U_r(ortho) U   = U_r
!   
!   as lu decomposition can be instable for orthogonalizing a set of 
!   vectors, the orthogonalization is refined iteratively. 
!   a conservative estimate of the orthogonalization error is used to
!   assess convergency
! 
!   a logical flag is set to false, if the routine fails, so that the
!   calling program can call a more robust orthogonalization routine
!   without aborting.
! 
!   arguments:
!   ==========
! 
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(inout) :: u_r, u_l
    logical,                   intent(inout) :: ok
! 
!   local variables:
!   ================
! 
    integer             :: it, it_micro
    real(dp)            :: error, dnrm2, alpha, unorm, shift
    logical             :: macro_done, micro_done, direct_solv
    real(dp), parameter :: tol_ortho_lu = 1.d-14 !two * epsilon(one)
    integer,  parameter :: maxit = 10
!
!   local scratch:
!   ==============
!
    real(dp), allocatable :: metric(:,:), msave(:,:), identity(:,:), ipiv(:)
! 
!   external functions:
!   ===================
!     
    external                                 :: dgemm, dgetrf, dtrtri
! 
!   get memory for the lu routine 
!
    allocate (metric(m,m), msave(m,m), identity(m,m))
    allocate (ipiv(m))
!
    metric      =  zero
    macro_done  = .false.
    direct_solv = .true.
!
    identity = zero
    forall(it = 1:m) identity(it, it) = one
!
!   start iteration and handle failure
!
    it = 0
    do while(.not. macro_done)
      it = it + 1
      if (it .gt. maxit) then
!
!       ortho_lu failed. return with an error message
!
        ok = .false.
        write(6,*) " maximum number of iterations reached in ortho_lu."
        write(6,'(A, E14.5)') " error of lu:", error
        return
      end if
!
!     assemble the metric 
!
      call dgemm('t','n',m,m,n,one,u_l,n,u_r,n,zero,metric,m)
      msave = metric
!
      print*, "iteration lu", it
      print *
      print *, "overlap of non biorthogonalized"
      call printMatrix(m,m,metric,m)
      print*
!
!
!     compute the lu factorization of the metric
!
      call dgetrf(m,m,metric,m,ipiv,info) 
!
!
!     if dgetrf failed, try a second time, after level-shifting the diagonal of the metric.
!
      if (info.ne.0) then
!  
        alpha      = 100.0_dp
        unorm      = dnrm2(n*m,u_r,1) ! TODO also do for u_l
        it_micro   = 0
        micro_done = .false.
! 
!       add larger and larger shifts to the diagonal until dgetrf manages to factorize it.
!  
        do while (.not. micro_done)
          it_micro = it_micro + 1
          if (it_micro.gt.maxit) then 
!  
!            ortho_lu failed. return with an error message
!  
            ok = .false.
            write(6,*) ' maximum number of the iterations reached.'
            stop
            return
          end if 
!  
          shift = max(epsilon(one)*alpha*unorm,tol_ortho)
          print*, "shift", shift
          metric = msave
          call diag_shift(m,shift,metric)
          print *
          print*, "Level shifted:"
          call printMatrix(m,m,metric,m)
          call dgetrf(m,m,metric,m,ipiv,info) 
          print *, "result"
          call printMatrix(m,m,metric,m)
          print*, "info", info
          print *
          alpha = alpha * 10.0_dp
          micro_done = info.eq.0
        end do
!  
      end if
!
      if (direct_solv) then
!
!       solve the matrix equations: U_l = U_l(ortho) l^t     
!                                   U_r = U_r(ortho) r
!
        call dtrsm('r','u','n','n',n,m,one,metric,m,u_r,n) 
        call dtrsm('r','l','t','u',n,m,one,metric,m,u_l,n) 
!
      else 
!
!       compute l^-1 and u^-1, using msave to store the inverse considering that
!       l is a unit triangular
!
        call dtrtri('u','n',m,metric,m,info)
        call checkInfo(info, "inverse of metric")
        call dtrtri('l','u',m,metric,m,info)
        call checkInfo(info, "inverse of metric")
!
        print *
        print *, "inverse of lower and upper and metric"
        call dgetrf(m,m,metric,m,ipiv,info) 
        print *
!
!       orthogonalize U_l and U_r:  U_l(ortho) = U_l l^(-t) 
!                                   U_r(ortho) = U_r u^(-1)  
!
        call dtrmm('r','u','n','n',n,m,one,metric,m,u_r,n) 
        call dtrmm('r','l','t','u',n,m,one,metric,m,u_l,n) 
!
      end if
!
!     calculate overlap of L^T & R, substract identity and check it's norm to
!     obtain the error
!
      call dgemm('t','n',m,m,n,one,u_l,n,u_r,n,zero,metric,m)
!
      metric = metric - identity
      error = dnrm2(m,metric,1)
      print*
      print *, "error", error, error .lt. tol_ortho_lu
      macro_done = error .lt. tol_ortho_lu
!
      print*, "print L^T*R overlap - I (residual overlap)"
      call printMatrix(m,m,metric,m)
      print*
!
    end do 
!
    ok = .true.
!
    deallocate (ipiv)
    deallocate (metric)
    deallocate (msave)
!
    return
  end subroutine ortho_lu
!
  subroutine biortho_vs_x(n,m,k,xl,xr,ul,ur)
    implicit none
    integer,                  intent(in)    :: n, m, k
    real(dp), dimension(n,m), intent(in)    :: xl, xr
    real(dp), dimension(n,k), intent(inout) :: ul, ur
!
!   local variables:
!
    integer               :: it, istat
    real(dp)              :: xu_norm(2), xx
    logical               :: done, ok
    real(dp), allocatable :: xu(:,:)
!
    integer, parameter    :: maxit = 20
    logical, parameter    :: allsvd = .true.
    real(dp)              :: dnrm2
!
    allocate (xu(m,k), stat = istat)
    call check_mem(istat)
!
    done = .false.
    it   = 0
!
    do while(.not. done)
      it = it + 1
      if (it.gt.maxit) stop 'biortho_vs_x failed.'
!
!     biorthogonalize ul and ur to xr and xl:
!
      call dgemm('t','n',m,k,n,one,xl,n,ur,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,xr,n,xu,m,one,ur,n)
      call dgemm('t','n',m,k,n,one,xr,n,ul,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,xl,n,xu,m,one,ul,n)
!
!     now, orthogonalize ur and ul. 
!     if allsvd is true, biorthogonalize them (to be tested)
!
      if (allsvd) then 
        call svd_biortho(n,k,ul,ur)
      else
        call ortho_cd(n,k,ul,xx,ok)
        call ortho_cd(n,k,ur,xx,ok)
      end if
!
!     compute the overlap between the orthonormalized ul, ur and xr, xl
!
      call dgemm('t','n',m,k,n,one,xl,n,ur,n,zero,xu,m)
      xu_norm(1) = dnrm2(m*k,xu,1)
      call dgemm('t','n',m,k,n,one,xr,n,ul,n,zero,xu,m)
      xu_norm(2) = dnrm2(m*k,xu,1)
!
      done = xu_norm(1) .lt. tol_ortho .and. xu_norm(2) .lt. tol_ortho
    end do
!
!   if we didn't make the vectors already biorthogonal, do it now:
!
    if (.not. allsvd) call svd_biortho(n,k,ul,ur)
!
    deallocate (xu)
    return
  end subroutine biortho_vs_x
!
  subroutine svd_biortho(n,m,u_l,u_r)
    implicit none
!
!   given two set of vectors, biorthogonalize them by computing the LU decomposition
!   of the overlap matrix and then solving
!
!     u_l = u_l l^t
!     u_r = u_r u
!
    integer,                  intent(in)    :: n, m
    real(dp), dimension(n,m), intent(inout) :: u_l, u_r
!
    integer               :: i, istat
    real(dp)              :: fac
!
    real(dp), allocatable :: over(:,:), u(:,:), s(:), vt(:,:), tmp(:,:)
!
    real(dp)              :: dnrm2
    external              :: dnrm2
!
!   allocate memory.
!
    allocate (over(m,m), s(m), u(m,m), vt(m,m), tmp(n,m), stat = istat)
    call check_mem(istat)
!
!   compute the overlap:
!
    call dgemm('t','n',m,m,n,one,u_l,n,u_r,n,zero,over,m)
!
!   compute its singular value decomposition:
!
    call dgesvd('a','a',m,m,over,m,s,u,m,vt,m,work,lwork,info)
!
!   compute l*u and r*v
!
    call dgemm('n','n',n,m,m,one,u_l,n,u,m,zero,tmp,n)
    u_l = tmp
    call dgemm('n','t',n,m,m,one,u_r,n,vt,m,zero,tmp,n)
    u_r = tmp
!
!   scale with square root of singular values
!   here, dropping redundant vectors could be a good idea...
!
    do i = 1, m
      fac = one/sqrt(s(i))
      u_l(:,i) = fac * u_l(:,i) 
      u_r(:,i) = fac * u_r(:,i) 
    end do
    deallocate(over, u, s, vt, tmp, stat = istat)
    return
  end subroutine svd_biortho
!
  real(dp) function norm_est(m,a)
!
!   compute a cheap estimate of the norm of a lower triangular matrix.
!   let a = d + o, where d = diag(a). as
!   
!   || a || <= || d || + || o ||
!
!   we compute || d || as max_i |d(i)| and || o || as its frobenius norm.
!   this is tight enough, and goes to 1 when a approaches the identity.
!
    implicit none
    integer,                  intent(in) :: m
    real(dp), dimension(m,m), intent(in) :: a
!
    integer  :: i, j
    real(dp) :: diag_norm, od_norm
!
    diag_norm = zero
    do i = 1, m
      diag_norm = max(diag_norm,abs(a(i,i)))
    end do
!
    od_norm = zero
    do i = 1, m
      do j = 1, i - 1
        od_norm = od_norm + a(i,j)**2
      end do
    end do
    od_norm = sqrt(od_norm)
!
    norm_est = diag_norm + od_norm
    return
  end function norm_est
!
  subroutine ortho_vs_x(n,m,k,x,u,ax,au)
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
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x, ax
    real(dp),  dimension(n,k), intent(inout) :: u, au
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it
    real(dp)               :: xu_norm, growth
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
    logical, parameter     :: useqr = .false.
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
    if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
    if (.not. ok .or. useqr) call ortho(n,k,u,au)
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
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
      if (.not. ok .or. useqr) call ortho(n,k,u,au)
!
!     the orthogonalization has introduced an error that makes the new
!     vector no longer fully orthogonal to x. assuming that u was 
!     orthogonal to x to machine precision before, we estimate the 
!     error with growth * eps, where growth is the product of the norms
!     of all the linear transformations applied to u.
!     if ortho_cd has failed, we just compute the overlap and its norm.
!
      if (.not. ok .or. useqr) then
        call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
        xu_norm = dnrm2(m*k,xu,1)
      else
        xu_norm = growth * epsilon(one)
      end if
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
    real(dp)               :: xu_norm, growth, xx(1)
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
    logical, parameter     :: useqr = .false.
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
    if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
    if (.not. ok .or. useqr) call ortho(n,k,u,xx)
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
      if (.not. useqr) call ortho_cd(n,k,u,growth,ok)
      if (.not. ok .or. useqr) call ortho(n,k,u,xx)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
!     note that, if we use ortho_cd, we estimate the norm of the overlap
!     using the growth factor returned in growth. 
!     see ortho_vs_x for more information.
!
      if (.not. ok .or. useqr) then
        call dgemm('t','n',m,k,n,one,bx,n,u,n,zero,xu,m)
        xu_norm = dnrm2(m*k,xu,1)
      else
        xu_norm = growth * epsilon(one)
      end if
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
    call ortho_vs_x(len_u,n_max,n_act,u_x,u_p,xx,xx)
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
    real(dp)              :: fac, diag_norm, out_norm, growth, xx(1)
    logical               :: ok
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
      call ortho_cd(n,m,evec,growth,ok)
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
        call ortho_cd(n,m,evec,growth,ok)
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
      if (iform.eq.1) write(6,100) mat(i,1:n)
      if (iform.eq.2) write(6,200) mat(i,1:n)
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
