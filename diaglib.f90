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
  private
!
! useful constants
!
  real(dp),    parameter   :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
  complex(dp), parameter   :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp), ctwo = (2.0_dp,0.0_dp) 
!
! convergence thresholds for orthogonalization
!
  real(dp), parameter    :: tol_ortho = 1.0e-12_dp
  !real(dp), parameter    :: tol_ortho = 1.0e2_dp
! 
! memory and info for lapack routines
!
  integer                   :: lwork, info
  real(dp),    allocatable  :: work(:), tau(:), rwork(:)
  complex(dp), allocatable  :: cwork(:), ctau(:)
!
! timings:
!
  real(dp)               :: t1(2), t2(2), t_diag(2), t_ortho(2), &
                            t_mv(2), t_tot(2)
!
! subroutines:
! ============
!
  public :: lobpcg_driver, davidson_driver, davidson_complex_driver, davidson_newcomplex_driver, gen_david_driver, & 
            caslr_driver, caslr_complex_driver, caslr_eff_driver, caslr_complex_eff_driver, caslr_newcomplex_driver, & 
            caslr_newcomplex_eff_driver, caslr_2x2_driver, ortho, b_ortho, ortho_cd, ortho_vs_x, ortho_vs_x_complex, &
            b_ortho_vs_x, b_ortho_vs_x_complex, b_ortho_complex, ortho_cd_complex, diag_shift_complex, & 
            ortho_complex, prtmat, prtmat_complex, check_sym_4n, ortho_vs_x_4n, ortho_cd_newcomplex
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
    real(dp)              :: sqrtn, tol_rms, tol_max
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
    call ortho_cd(.false.,n,n_max,vp,xx,ok)
    call ortho_cd(.false.,n,n_max,vm,xx,ok)
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
!fl
!         write(6,*) 'eig: ', eig(i_eig)
!         write(6,*) 'up ', i_eig
!         write(6,'(10f12.6)') up(1:ldu,i_eig)
!         write(6,*) 'um ', i_eig
!         write(6,'(10f12.6)') um(1:ldu,i_eig)
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
!fl
!         write(6,*) 'eig: ', eig(i_eig)
!         write(6,*) 'up ', i_eig
!         write(6,'(10f12.6)') up(1:ldu,i_eig)
!         write(6,*) 'um ', i_eig
!         write(6,'(10f12.6)') um(1:ldu,i_eig)
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
        call ortho_vs_x(.false.,n,ldu,n_act,vp,vp(1,i_beg),xx,xx)
        call ortho_vs_x(.false.,n,ldu,n_act,vm,vm(1,i_beg),xx,xx)
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
        call ortho_cd(.false.,n,n_max,vp,xx,ok)
        call ortho_cd(.false.,n,n_max,vm,xx,ok)
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
  subroutine caslr_2x2_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                              lambdamul,omegamul,lrprec,eig,evec,ok)
!
!   fake olsen driver: b(+) = (b+ b+), b(-) = (b- -b-)
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
    external                                     :: lambdamul, omegamul, &
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
    real(dp)              :: sqrtn, tol_rms, tol_max
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
    real(dp), allocatable :: vp2(:,:), vm2(:,:), lvp2(:,:), lvm2(:,:), bvp2(:,:), bvm2(:,:)
    real(dp), allocatable :: rp2(:,:), rm2(:,:), rr2(:,:), temp(:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp), allocatable :: up(:,:), um(:,:), eigp2(:,:), eigm2(:,:), bp2(:,:), bm2(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: epmat(:,:), emmat(:,:), smat(:,:)
! 
!   test full 4x4 problem
!
!    real(dp), allocatable :: bp_r(:,:), bp_i(:,:), bm_r(:,:), bm_i(:,:), evec_old(:,:)
!    real(dp), allocatable :: overlap(:,:), vp_old(:,:), overlap_vp(:,:)
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
    allocate (vp(n/2,lda), vm(n/2,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), &
              vp2(n2,lda), vm2(n2,lda), lvp2(n2,lda), lvm2(n2,lda), bvp2(n2,lda), bvm2(n2,lda), &
              rp2(n2,n_max), rm2(n2,n_max), rr2(n2,n_max), temp(n2), stat = istat)
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
    allocate (up(lda,n_max), um(lda,n_max), eigp2(n2,n_max), eigm2(n2,n_max), &
              bp2(n2,n_max), bm2(n2,n_max), stat = istat)
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n2,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag     = zero
    t_ortho    = zero
    t_mv       = zero
    t_tot      = zero
    vp         = zero
    vm         = zero
    bvp        = zero
    bvm        = zero
    lvp        = zero
    lvm        = zero
    vp2        = zero
    vm2        = zero
    bvp2       = zero
    bvm2       = zero
    lvp2       = zero
    lvm2       = zero
    a_red      = zero
    ok         = .false.
    done       = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
!   keep in mind that: X_eig   = (Yr  Zr  Yi -Zi)
!
    do i_eig = 1, n_max, 2
      !eig i
      vp(:,i_eig)   = (evec(1:n/2,i_eig) + evec(n/2+1:n,i_eig))/2.0_dp
      vm(:,i_eig)   = (evec(1:n/2,i_eig) - evec(n/2+1:n,i_eig))/2.0_dp
      !eig i+1
      vp(:,i_eig+1) = (+evec(n+1:(3*n/2),i_eig) + evec((3*n/2)+1:n2,i_eig))/2.0_dp
      vm(:,i_eig+1) = (+evec(n+1:(3*n/2),i_eig) - evec((3*n/2)+1:n2,i_eig))/2.0_dp
    end do
!
    vp2(1:n/2,:)         = vp
    vp2(n/2+1:n,:)       = vp
    vm2(1:n/2,:)         = vm
    vm2(n/2+1:n,:)       = -vm
    !
    vp = zero
    vm = zero
    !
    do i_eig = 1, n_max, 2
      !eig
      vp(:,i_eig)   = (evec(n+1:(3*n)/2,i_eig) - evec((3*n/2)+1:n2,i_eig))/2.0_dp
      vm(:,i_eig)   = (evec(n+1:(3*n)/2,i_eig) + evec((3*n/2)+1:n2,i_eig))/2.0_dp
      !eig i+1
      vp(:,i_eig+1) = (-evec(1:n/2,i_eig) + evec(n/2+1:n,i_eig))/2.0_dp
      vm(:,i_eig+1) = (-evec(1:n/2,i_eig) - evec(n/2+1:n,i_eig))/2.0_dp
    end do
    !
    vp2(n+1:(3*n2)/4,:)  = vp
    vp2((3*n2)/4+1:n2,:) = -vp
    vm2(n+1:(3*n2)/4,:)  = vm
    vm2((3*n2)/4+1:n2,:) = vm
!
    !call check_sym_4n(1,n2,n_max,vp2)
    !call check_sym_4n(2,n2,n_max,vm2)
    call ortho_cd(.false.,n2,n_max,vp2,xx,ok)
    call ortho_cd(.false.,n2,n_max,vm2,xx,ok)
    !call check_sym_4n(1,n2,n_max,vp2)
    !call check_sym_4n(2,n2,n_max,vm2)
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
      call lambdamul(n2,n_act,vp2(1,i_beg),lvp2(1,i_beg))
      call lambdamul(n2,n_act,vm2(1,i_beg),lvm2(1,i_beg))
      call omegamul(n2,n_act,vp2(1,i_beg),bvm2(1,i_beg))
      call omegamul(n2,n_act,vm2(1,i_beg),bvp2(1,i_beg))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,ldu,n2,one,vp2,n2,lvp2,n2,zero,epmat,lda)
      call dgemm('t','n',ldu,ldu,n2,one,vm2,n2,lvm2,n2,zero,emmat,lda)
      call dgemm('t','n',ldu,ldu,n2,one,vm2,n2,bvm2,n2,zero,smat,lda)
!
      e_red = zero
      a_red = zero
      s_red = zero
      a_red(1:ldu,1:ldu)             = epmat(1:ldu,1:ldu)
      a_red(ldu+1:2*ldu,ldu+1:2*ldu) = emmat(1:ldu,1:ldu)
      s_red(1:ldu,ldu+1:2*ldu)       = transpose(smat(1:ldu,1:ldu))
      s_red(ldu+1:2*ldu,1:ldu)       = smat(1:ldu,1:ldu)
!
      call get_time(t1)
!
!     default algorithm: solve the 2n-dimensional inverse problem.
!
      call dsygv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,work,lwork,info)
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      do i_eig = 1, n_max
        eig(i_eig)      = one/e_red(2*ldu - i_eig + 1)
        up(1:ldu,i_eig) = s_red(1:ldu,2*ldu - i_eig + 1)
        um(1:ldu,i_eig) = s_red(ldu+1:2*ldu,2*ldu - i_eig + 1)
      end do
!
      call get_time(t2)
      t_diag = t_diag + t2 - t1
      call dgemm('n','n',n2,n_max,ldu,one,vp2,n2,up,lda,zero,eigp2,n2)
      call dgemm('n','n',n2,n_max,ldu,one,vm2,n2,um,lda,zero,eigm2,n2)
!
      do i_eig = 1, n_max
        evec(:,i_eig) = eigp2(:,i_eig) + eigm2(:,i_eig)
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n2,n_max,ldu,one,lvp2,n2,up,lda,zero,rp2,n2)
      call dgemm('n','n',n2,n_max,ldu,one,lvm2,n2,um,lda,zero,rm2,n2)
      call dgemm('n','n',n2,n_max,ldu,one,bvp2,n2,um,lda,zero,bp2,n2)
      call dgemm('n','n',n2,n_max,ldu,one,bvm2,n2,up,lda,zero,bm2,n2)
!      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        call daxpy(n2,-eig(i_eig),bp2(:,i_eig),1,rp2(:,i_eig),1)
        call daxpy(n2,-eig(i_eig),bm2(:,i_eig),1,rm2(:,i_eig),1)
        r_norm(1,i_eig) = sqrt(dot_product(rp2(:,i_eig),rp2(:,i_eig)) + &
                               dot_product(rm2(:,i_eig),rm2(:,i_eig)))/sqrtn
        r_norm(2,i_eig) = maxval(abs(rp2(:,i_eig) + rm2(:,i_eig)))
      end do
!
!     check the sign of the first element to avoid broken-symmetry failure
!
!      temp = zero
!      do i_eig = 1, n_targ, 2
!        if ((rp2(1,i_eig)*rm2(n+1,i_eig+1)) .lt. zero) then
!          write(6,*) 'I swapped what had to be swapped'
!          temp           = rp2(:,i_eig)
!          rp2(:,i_eig)   = rm2(:,i_eig+1)
!          rm2(:,i_eig+1) = temp
!        end if 
!        temp             = zero
!        if ((rm2(1,i_eig)*rp2(n+1,i_eig+1)) .lt. zero) then
!          temp           = rm2(:,i_eig)
!          rm2(:,i_eig)   = rp2(:,i_eig+1)
!          rp2(:,i_eig+1) = temp
!        end if 
!      end do
!      
!      do i_eig = 1, n_targ, 2
!        write(6,'(4(a3,i1,11x))') 'R+_', i_eig, 'R+_', i_eig+1, 'R-_', i_eig, 'R-_', i_eig+1
!        do i = 1, 10
!          !write(6,'(4f15.8)') rp2(i,i_eig), rp2(n+i,i_eig+1), rm2(i,i_eig), rm2(n+i,i_eig+1)
!          write(6,'(2f15.8)') rp2(i,i_eig)+rm2(i,i_eig), rp2(n+i,i_eig+1)+rm2(n+i,i_eig+1)
!        end do
!      end do
!
!     check convergence. lock the first contiguous converged eigenvalues
!     by setting the logical array "done" to true.
!
      do i_eig = 1, n_targ, 2
        if (done(i_eig)) cycle
        done(i_eig)     = r_norm(1,i_eig).lt.tol_rms .and. &
                          r_norm(2,i_eig).lt.tol_max .and. &
                          it.gt.1
        done(i_eig+1)   = r_norm(1,i_eig+1).lt.tol_rms .and. &
                          r_norm(2,i_eig+1).lt.tol_max .and. &
                          it.gt.1
!
!     check coupled-eigenvalues convergece: if one eigval from the same couple
!     is false, both need to be put false (even if the other was .true.)
!
        if (done(i_eig) .neqv. done(i_eig+1)) then 
          done(i_eig)   = .false.
          done(i_eig+1) = .false.
        end if
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
        call lrprec(n2,n_act,eig(ind),rp2(1,ind),rm2(1,ind),vp2(1,i_beg),vm2(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call ortho_vs_x_4n(.false.,n2,ldu,n_act,vp2,vp2(1,i_beg),xx,xx,1)
        call ortho_vs_x_4n(.false.,n2,ldu,n_act,vm2,vm2(1,i_beg),xx,xx,2)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        restart = .true.
!
!       initialize indices back to their starting values 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!
        n_act = n_max 
        vp    = zero
        vm    = zero
        vp2   = zero 
        vm2   = zero 
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max, 2
          !eig i
          vp(:,i_eig)    = (evec(1:n/2,i_eig) + evec(n/2+1:n,i_eig))/2.0_dp
          vm(:,i_eig)    = (evec(1:n/2,i_eig) - evec(n/2+1:n,i_eig))/2.0_dp
          !eig i+1
          vp(:,i_eig+1) = (+evec(n+1:(3*n/2),i_eig) + evec((3*n/2)+1:n2,i_eig))/2.0_dp
          vm(:,i_eig+1) = (+evec(n+1:(3*n/2),i_eig) - evec((3*n/2)+1:n2,i_eig))/2.0_dp
        end do
!
!       build b(+) and b(-)
!    
        vp2(1:n/2,:)         =  vp
        vp2(n/2+1:n,:)       =  vp
        vm2(1:n/2,:)         =  vm
        vm2(n/2+1:n,:)       = -vm
        !
        vp = zero
        vm = zero
        !
        do i_eig = 1, n_max, 2
          !eig
          vp(:,i_eig)   = (evec(n+1:(3*n)/2,i_eig) - evec((3*n/2)+1:n2,i_eig))/2.0_dp
          vm(:,i_eig)   = (evec(n+1:(3*n)/2,i_eig) + evec((3*n/2)+1:n2,i_eig))/2.0_dp
          !eig i+1
          vp(:,i_eig+1) = (-evec(1:n/2,i_eig) + evec(n/2+1:n,i_eig))/2.0_dp
          vm(:,i_eig+1) = (-evec(1:n/2,i_eig) - evec(n/2+1:n,i_eig))/2.0_dp
        end do
        !
        vp2(n+1:(3*n2)/4,:)  =  vp
        vp2((3*n2)/4+1:n2,:) = -vp
        vm2(n+1:(3*n2)/4,:)  =  vm
        vm2((3*n2)/4+1:n2,:) =  vm
!        
!test        do i_eig = 1, n_max, 2
!test          write(6,*) 'evec-vp-vm; i_eig: ', i_eig
!test          do i = 1, n2
!test            if (i == n/2+1)   write(6,*) '50'
!test            if (i == n+1)     write(6,*) '100'
!test            if (i == 3*n/2+1) write(6,*) '150'
!test            write(6,'(2f20.16)') evec(i,i_eig) - vp2(i,i_eig) - vm2(i,i_eig), evec(i,i_eig+1) - vp2(i,i_eig+1) - vm2(i,i_eig+1)
!test            if (i == n2)      write(6,*) '200'
!test          end do
!test        end do
!    
        call ortho_cd(.false.,n2,n_max,vp2,xx,ok)
        call ortho_cd(.false.,n2,n_max,vm2,xx,ok)
!
        lvp   = zero
        lvm   = zero
        bvp   = zero
        bvm   = zero
        lvp2  = zero
        lvm2  = zero
        bvp2  = zero
        bvm2  = zero 
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
               vp2,vm2,lvp2,lvm2,bvp2,bvm2,rp2,rm2,rr2, &
               epmat,emmat,smat,up,um,eigp2,eigm2,bp2,bm2)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_2x2_driver
!
  subroutine caslr_newcomplex_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                                     apbmul,ambmul,spdmul,smdmul,lrprec,eig,evecre,evecim,ok)
    use utils
    implicit none
    logical, intent(in)                          :: verbose
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n_max),    intent(inout) :: eig
    real(dp), dimension(n2,n_max), intent(inout) :: evecre, evecim
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
    real(dp)              :: sqrtn, sqrt2n, tol_rms, tol_max
    real(dp)              :: full_norm
    real(dp)              :: xx(1), lw_svd(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
!   real
    real(dp), allocatable :: vpre(:,:), vmre(:,:), lvpre(:,:), lvmre(:,:), bvpre(:,:), bvmre(:,:)
    real(dp), allocatable :: rpre(:,:), rmre(:,:), rrre(:,:), r_norm(:,:)
!   imaginary
    real(dp), allocatable :: vpim(:,:), vmim(:,:), lvpim(:,:), lvmim(:,:), bvpim(:,:), bvmim(:,:)
    real(dp), allocatable :: rpim(:,:), rmim(:,:), rrim(:,:)
!
!   debug
!
    real(dp), allocatable :: vpp(:,:), vmm(:,:) 
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
!   real
    real(dp), allocatable :: up(:,:), um(:,:), eigpre(:,:), eigmre(:,:), bpre(:,:), bmre(:,:)
    real(dp), allocatable :: eigpre_(:,:), eigmre_(:,:)
!   imaginary
    real(dp), allocatable :: eigpim(:,:), eigmim(:,:), bpim(:,:), bmim(:,:)
    real(dp), allocatable :: eigpim_(:,:), eigmim_(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: epmat(:,:),   emmat(:,:),   smat(:,:)
    real(dp), allocatable :: epmatre(:,:), emmatre(:,:), smatre(:,:)
    real(dp), allocatable :: epmatim(:,:), emmatim(:,:), smatim(:,:)
!    
    real(dp), allocatable :: epmmatre(:,:), epmmatim(:,:), epmmat(:,:)
    real(dp), allocatable :: empmatre(:,:), empmatim(:,:), empmat(:,:) 
    real(dp), allocatable :: spmatre(:,:), spmatim(:,:), spmat(:,:)   
    real(dp), allocatable :: smmatre(:,:), smmatim(:,:), smmat(:,:)   
!
!   reduced 4nx4n problem matrices: 
!   p stands for ++, m stands for --, r stands for real, i stand for imaginary
!
!    real(dp), allocatable :: ep_rr(:,:), ep_ri(:,:), ep_ir(:,:), ep_ii(:,:)
!    real(dp), allocatable :: em_rr(:,:), em_ri(:,:), em_ir(:,:), em_ii(:,:)
!    real(dp), allocatable :: smp_rr(:,:), smp_ri(:,:), smp_ir(:,:), smp_ii(:,:)
!    real(dp), allocatable :: spm_rr(:,:), spm_ri(:,:), spm_ir(:,:), spm_ii(:,:)
!    real(dp), allocatable :: emat4(:,:), smat4(:,:), eigen4(:) 
    
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
    !dim_dav = 5
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
    allocate (vpre(n,lda), vpim(n,lda), vmre(n,lda), vmim(n,lda), lvpre(n,lda), lvpim(n,lda), lvmre(n,lda), lvmim(n,lda), &
              bvpre(n,lda), bvpim(n,lda), bvmre(n,lda), bvmim(n,lda), &
              rpre(n,n_max), rpim(n,n_max), rmre(n,n_max), rmim(n,n_max), rrre(n,n_max), rrim(n,n_max), stat = istat)
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
              e_red(lda2), epmat(lda,lda), epmatre(lda,lda), epmatim(lda,lda), emmat(lda,lda), &
              emmatre(lda,lda), emmatim(lda,lda), smat(lda,lda), smatre(lda,lda), smatim(lda,lda), & 
!      
              epmmatre(lda,lda), epmmatim(lda,lda), epmmat(lda,lda), & 
              empmatre(lda,lda), empmatim(lda,lda), empmat(lda,lda), &
              spmatre(lda,lda), spmatim(lda,lda), spmat(lda,lda), &
              smmatre(lda,lda), smmatim(lda,lda), smmat(lda,lda), & 
!
!     reduced 4nx4n problem allocation
!
!              ep_rr(lda,lda), ep_ri(lda,lda), ep_ir(lda,lda), ep_ii(lda,lda), &
!              em_rr(lda,lda), em_ri(lda,lda), em_ir(lda,lda), em_ii(lda,lda), &
!              spm_rr(lda,lda), spm_ri(lda,lda), spm_ir(lda,lda), spm_ii(lda,lda), &
!              smp_rr(lda,lda), smp_ri(lda,lda), smp_ir(lda,lda), smp_ii(lda,lda), & 
!              emat4(4*lda,4*lda), smat4(4*lda,4*lda), eigen4(4*lda), & 
               stat = istat) 
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigpre(n,n_max), eigpim(n,n_max), eigmre(n,n_max), eigmim(n,n_max), & 
              bpre(n,n_max), bpim(n,n_max), bmre(n,n_max), bmim(n,n_max), & 
              eigpre_(n,n_max), eigmre_(n,n_max), eigpim_(n,n_max), eigmim_(n,n_max),stat = istat)
    call check_mem(istat)
!
!   debug
!
    allocate (vpp(2*n,lda), vmm(2*n,lda), stat = istat)
    call check_mem(istat)
!
    write(6,*) 'EFFICIENT COMPLEX OLSEN IMPLEMENTATION'
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn    = sqrt(real(n,dp))
    sqrt2n   = sqrt(real(n2,dp))
    tol_rms  = tol
    tol_max  = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag    = zero
    t_ortho   = zero
    t_mv      = zero
    t_tot     = zero
    vpre      = zero
    vpim      = zero
    vmre      = zero
    vmim      = zero
    bvpre     = zero
    bvpim     = zero
    bvmre     = zero
    bvmim     = zero
    lvpre     = zero
    lvpim     = zero
    lvmre     = zero
    lvmim     = zero
!    emat4     = zero
!    smat4     = zero
    a_red     = zero
    ok        = .false.
    done      = .false.
!
!   debug
!
    vpp = zero
    vmm = zero
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
!
!   keep in mind that: X_eig   = (Yr  Zr  Yi -Zi)
!                      X_eig+1 = (-Yi Zi Yr Zr)
!                      evecre  = (Yr Zr), (-Yi Zi)
!                      evecim  = (Yi -Zi), (Yr Zr) 
!
    do i_eig = 1, n_max, 2
!     real 
      vpre(:,i_eig)   = (evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
      vmre(:,i_eig)   = (evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
      vpre(:,i_eig+1) = (+evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vmre(:,i_eig+1) = (+evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!     imaginary 
      vpim(:,i_eig)   = (evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig)   = (evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vpim(:,i_eig+1) = (-evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig+1) = (-evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
    end do
!
    vpp(1:n,:)     = vpre
    vpp(n+1:n2,:)  = vpim
    vmm(1:n,:)     = vmre
    vmm(n+1:n2,:)  = vmim
    call ortho_cd_newcomplex(.false.,n2,n_max,vpp,xx,ok)
    call ortho_cd_newcomplex(.false.,n2,n_max,vmm,xx,ok)
    vpre = vpp(1:n,:)
    vpim = vpp(n+1:n2,:)
    vmre = vmm(1:n,:)
    vmim = vmm(n+1:n2,:)
!   real 
!    call ortho_cd(.false.,n,n_max,vpre,xx,ok)
!    call ortho_cd(.false.,n,n_max,vmre,xx,ok)
!   imaginary
!    call ortho_cd(.false.,n,n_max,vpim,xx,ok)
!    call ortho_cd(.false.,n,n_max,vmim,xx,ok)
!    do i_eig = 1, n_max
!      write(6,*) 'vp 1 = r+ r+ i+ -i+', i_eig
!      do i = 1, n
!        write(6,'(4f8.4)') vpre(i,i_eig), vpre(i,i_eig), vpim(i,i_eig), -vpim(i,i_eig) 
!      end do
!      write(6,*) 'vm 2 = -i+ i+ r+ r+'
!      do i = 1, n
!        write(6,'(4f8.4)') vmre(i,i_eig), -vmre(i,i_eig), vmim(i,i_eig), vmim(i,i_eig) 
!      end do
!    end do
!
!    vpp(1:n,:)    = vpre 
!    vpp(n+1:n2,:) = vpim 
!    vmm(1:n,:)    = vmre
!    vmm(n+1:n2,:) = vmim
!    write(6,*) 'dotprod', dot_product(vpp(:,1),vpp(:,2))
!    write(6,*) 'dotprod', dot_product(vmm(:,1),vmm(:,2))
!    write(6,*) 'Vr+*Vi+', dot_product(vpre(:,1),vpim(:,1))
!    write(6,*) 'Vr-*Vi-', dot_product(vmre(:,1),vmim(:,1))
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
!
!     compute the matrix-vector products
!
      call apbmul(0,n,n_act,vpre(1,i_beg),vpim(1,i_beg),lvpre(1,i_beg))   !sigma_re+
      call apbmul(1,n,n_act,vpre(1,i_beg),vpim(1,i_beg),lvpim(1,i_beg))   !sigma_im+
      call ambmul(0,n,n_act,vmre(1,i_beg),vmim(1,i_beg),lvmre(1,i_beg))   !sigma_re-
      call ambmul(1,n,n_act,vmre(1,i_beg),vmim(1,i_beg),lvmim(1,i_beg))   !sigma_im-
      call spdmul(0,n,n_act,vpre(1,i_beg),vpim(1,i_beg),bvmre(1,i_beg))   !tau_re-
      call spdmul(1,n,n_act,vpre(1,i_beg),vpim(1,i_beg),bvmim(1,i_beg))   !tau_im-
      call smdmul(0,n,n_act,vmre(1,i_beg),vmim(1,i_beg),bvpre(1,i_beg))   !tau_re+
      call smdmul(1,n,n_act,vmre(1,i_beg),vmim(1,i_beg),bvpim(1,i_beg))   !tau_im+
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
!     reduced 4nx4n problem implementation 
!
!new      call dgemm('t','n',ldu,ldu,n,one,vpre,n,lvpre,n,zero,ep_rr,lda)  !Err_++=(Vre+)^T*sigma_re+
!new!new      write(6,*) 'Err++'
!new!new      print '(4f10.4)', transpose(ep_rr(1:ldu,1:ldu))
!new      call dgemm('t','n',ldu,ldu,n,one,vpre,n,lvpim,n,zero,ep_ri,lda)  !Eri_++=(Vre+)^T*sigma_im+
!new      call dgemm('t','n',ldu,ldu,n,one,vpim,n,lvpre,n,zero,ep_ir,lda)  !Eir_++=(Vim+)^T*sigma_re+
!new      call dgemm('t','n',ldu,ldu,n,one,vpim,n,lvpim,n,zero,ep_ii,lda)  !Eii_++=(Vim+)^T*sigma_im+
!new      !
!new      call dgemm('t','n',ldu,ldu,n,one,vmre,n,lvmre,n,zero,em_rr,lda)  !Err_--=(Vre-)^T*sigma_re-
!new      call dgemm('t','n',ldu,ldu,n,one,vmre,n,lvmim,n,zero,em_ri,lda)  !Eri_--=(Vre-)^T*sigma_im-
!new      call dgemm('t','n',ldu,ldu,n,one,vmim,n,lvmre,n,zero,em_ir,lda)  !Eir_--=(Vim-)^T*sigma_re-
!new      call dgemm('t','n',ldu,ldu,n,one,vmim,n,lvmim,n,zero,em_ii,lda)  !Eii_--=(Vim-)^T*sigma_im-
!new      !
!new      call dgemm('t','n',ldu,ldu,n,one,vpre,n,bvpre,n,zero,spm_rr,lda) !Srr_+-=(Vre+)^T*tau_re+
!new!new      print *, 'Srr+-'
!new!new      print '(4f10.4)', transpose(spm_rr(1:4,1:4))
!new      call dgemm('t','n',ldu,ldu,n,one,vpre,n,bvpim,n,zero,spm_ri,lda) !Sri_+-=(Vre+)^T*tau_im+
!new!new      print *, 'Sri+-'
!new!new      print '(4f10.4)', transpose(spm_ri(1:4,1:4))
!new      call dgemm('t','n',ldu,ldu,n,one,vpim,n,bvpre,n,zero,spm_ir,lda) !Sir_+-=(Vim+)^T*tau_re+
!new!new      print *, 'Sir+-'
!new!new      print '(4f10.4)', transpose(spm_ir(1:4,1:4))
!new      call dgemm('t','n',ldu,ldu,n,one,vpim,n,bvpim,n,zero,spm_ii,lda) !Sii_+-=(Vim+)^T*tau_im+
!new!new      print *, 'Sii+-'
!new!new      print '(4f10.4)', transpose(spm_ii(1:4,1:4))
!new      !
!new      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvmre,n,zero,smp_rr,lda) !Srr_-+=(Vre-)^T*tau_re-
!new!new      print *, 'Srr-+'
!new!new      print '(4f10.4)', transpose(smp_rr(1:4,1:4))
!new      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvmim,n,zero,smp_ri,lda) !Sri_-+=(Vre-)^T*tau_im-
!new!new      print *, 'Sri-+'
!new!new      print '(4f10.4)', transpose(smp_ri(1:4,1:4))
!new      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvmre,n,zero,smp_ir,lda) !Sir_-+=(Vim-)^T*tau_re-
!new!new      print *, 'Sir-+'
!new!new      print '(4f10.4)', transpose(smp_ir(1:4,1:4))
!new      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvmim,n,zero,smp_ii,lda) !Sii_-+=(Vim-)^T*tau_im-
!new!new      print *, 'Sii-+'
!new!new      print '(4f10.4)', transpose(smp_ii(1:4,1:4))
!new!      
!new!     now, build the 4nx4n reduced matrices E and S 
!new!
!new      emat4(1:ldu,1:ldu)                 = ep_rr(1:ldu,1:ldu) 
!new      emat4(1:ldu,ldu+1:2*ldu)           = ep_ri(1:ldu,1:ldu)
!new      emat4(ldu+1:2*ldu,1:ldu)           = ep_ir(1:ldu,1:ldu)     
!new      emat4(ldu+1:2*ldu,ldu+1:2*ldu)     = ep_ii(1:ldu,1:ldu)
!new      emat4(2*ldu+1:3*ldu,2*ldu+1:3*ldu) = em_rr(1:ldu,1:ldu)
!new      emat4(2*ldu+1:3*ldu,3*ldu+1:4*ldu) = em_ri(1:ldu,1:ldu)
!new      emat4(3*ldu+1:4*ldu,2*ldu+1:3*ldu) = em_ir(1:ldu,1:ldu)
!new      emat4(3*ldu+1:4*ldu,3*ldu+1:4*ldu) = em_ii(1:ldu,1:ldu)  
!new      ! 
!new      smat4(1:ldu,2*ldu+1:3*ldu)         = spm_rr(1:ldu,1:ldu)  
!new      smat4(1:ldu,3*ldu+1:4*ldu)         = spm_ri(1:ldu,1:ldu)
!new      smat4(ldu+1:2*ldu,2*ldu+1:3*ldu)   = spm_ir(1:ldu,1:ldu)
!new      smat4(ldu+1:2*ldu,3*ldu+1:4*ldu)   = spm_ii(1:ldu,1:ldu)
!new      smat4(2*ldu+1:3*ldu,1:ldu)         = smp_rr(1:ldu,1:ldu)
!new      smat4(2*ldu+1:3*ldu,ldu+1:2*ldu)   = smp_ri(1:ldu,1:ldu)
!new      smat4(3*ldu+1:4*ldu,1:ldu)         = smp_ir(1:ldu,1:ldu)
!new      smat4(3*ldu+1:4*ldu,ldu+1:2*ldu)   = smp_ii(1:ldu,1:ldu)
!new      print *, 'E red'
!new      print '(16f8.4)', transpose(emat4(1:16,1:16))
!new      print *, 'S red'
!new      print '(16f8.4)', transpose(smat4(1:16,1:16))
! 
!dc      call dgemm('t','n',ldu,ldu,n,one,vpre,n,lvmre,n,zero,epmmatre,lda)  !Ere+-=(Vre+)^T*sigma_re-
!dc      call dgemm('t','n',ldu,ldu,n,one,vpim,n,lvmim,n,zero,epmmatim,lda)  !Eim+-=(Vim+)^T*sigma_im-
!dc      epmmat = epmmatre + epmmatim
!dc      !write(6,*) 'E+-', epmmat
!dc      !
!dc      call dgemm('t','n',ldu,ldu,n,one,vmre,n,lvpre,n,zero,empmatre,lda)  !Ere-+=(Vre-)^T*sigma_re+
!dc      call dgemm('t','n',ldu,ldu,n,one,vmim,n,lvpim,n,zero,empmatim,lda)  !Eim-+=(Vim-)^T*sigma_im+
!dc      empmat = empmatre + empmatim 
!dc      !write(5,*) 'E-+', empmat
!dc      !
!dc      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvpre,n,zero,spmatre,lda)    !S+re=(Vre-)^T*tau_re+
!dc      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvpim,n,zero,spmatim,lda)    !S+im=(Vim-)^T*tau_im+
!dc      spmat = spmatre + spmatim
!dc      write(6,*) 'S++', spmat
!dc      !
!dc      call dgemm('t','n',ldu,ldu,n,one,vpre,n,bvmre,n,zero,smmatre,lda)    !S-re=(Vre+)^T*tau_re-
!dc      call dgemm('t','n',ldu,ldu,n,one,vpim,n,bvmim,n,zero,smmatim,lda)    !S-im=(Vim+)^T*tau_im-
!dc      smmat = smmatre + smmatim
!dc      write(6,*) 'S--', smmat
!dc      stop
!
      call dgemm('t','n',ldu,ldu,n,one,vpre,n,lvpre,n,zero,epmatre,lda)   !Ere+=(Vre+)^T*sigma_re+
      call dgemm('t','n',ldu,ldu,n,one,vpim,n,lvpim,n,zero,epmatim,lda)   !Eim+=(Vim+)^T*sigma_im+
      epmat = two*(epmatre + epmatim)
      call dgemm('t','n',ldu,ldu,n,one,vmre,n,lvmre,n,zero,emmatre,lda)   !Ere-=(Vre-)^T*sigma_re-
      call dgemm('t','n',ldu,ldu,n,one,vmim,n,lvmim,n,zero,emmatim,lda)   !Eim-=(Vim-)^T*sigma_im-
      emmat = two*(emmatre + emmatim)
      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvmre,n,zero,smatre,lda)    !Sre=(Vre-)^T*tau_re-
      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvmim,n,zero,smatim,lda)    !Sim=(Vim-)^T*tau_im-
      smat = two*(smatre + smatim)
!
      e_red = zero
      a_red = zero
      s_red = zero
      a_red(1:ldu,1:ldu)             = epmat(1:ldu,1:ldu)
      a_red(ldu+1:2*ldu,ldu+1:2*ldu) = emmat(1:ldu,1:ldu)
      s_red(1:ldu,ldu+1:2*ldu)       = transpose(smat(1:ldu,1:ldu))
      s_red(ldu+1:2*ldu,1:ldu)       = smat(1:ldu,1:ldu)
!
      call get_time(t1)
!
!     default algorithm: solve the 2n-dimensional inverse problem.
!
      call dsygv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,work,lwork,info)
      print *, 'info dsygv', info
!
!     Olsen algorithm: solve the whole inverse problem (4nx4n) 
!
!new      call dsygv(1,'v','l',4*ldu,smat4,4*lda,emat4,4*lda,eigen4,work,lwork,info)
!new      write(6,*) 'info dsygv', info
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      do i_eig = 1, n_max
        eig(i_eig)      = one/e_red(2*ldu - i_eig + 1)
        up(1:ldu,i_eig) = s_red(1:ldu,2*ldu - i_eig + 1)
        um(1:ldu,i_eig) = s_red(ldu+1:2*ldu,2*ldu - i_eig + 1)
!new        print *, one/eigen4(i_eig)
!new        eig(i_eig)      = one/eigen4(4*ldu - i_eig + 1)
        !print *, i_eig, eig(i_eig)
      end do
!
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     real part 
!
      call dgemm('n','n',n,n_max,ldu,one,vpre,n,up,lda,zero,eigpre,n)
      call dgemm('n','n',n,n_max,ldu,one,vmre,n,um,lda,zero,eigmre,n)
!
!     imaginary part
!
      call dgemm('n','n',n,n_max,ldu,one,vpim,n,up,lda,zero,eigpim,n)
      call dgemm('n','n',n,n_max,ldu,one,vmim,n,um,lda,zero,eigmim,n)
!
!     normalize 
!
!      do i_eig = 1, n_max
!        full_norm = sqrt(two*(dot_product(eigpre(:,i_eig),eigpre(:,i_eig)) + & 
!                              dot_product(eigpim(:,i_eig),eigpim(:,i_eig))))
!        eigpre(:,i_eig) = eigpre(:,i_eig)/full_norm
!        eigpim(:,i_eig) = eigpim(:,i_eig)/full_norm
!        full_norm = sqrt(two*(dot_product(eigmre(:,i_eig),eigmre(:,i_eig)) + & 
!                              dot_product(eigmim(:,i_eig),eigmim(:,i_eig))))
!        eigmre(:,i_eig) = eigmre(:,i_eig)/full_norm
!        eigmim(:,i_eig) = eigmim(:,i_eig)/full_norm
!      end do
!
      do i_eig = 1, n_max!, 2
!
!       being the evec defined as Xeig   = (Yr Zr Yi -Zi)
!                                 Xeig+1 = (-Yi Zi Yr Zr)
!
        evecre(1:n,i_eig)    =  eigpre(:,i_eig) + eigmre(:,i_eig)   !Yre
        evecre(n+1:n2,i_eig) =  eigpre(:,i_eig) - eigmre(:,i_eig)   !Zre
        evecim(1:n,i_eig)    =  eigpim(:,i_eig) + eigmim(:,i_eig)   !Yim
        evecim(n+1:n2,i_eig) = -eigpim(:,i_eig) + eigmim(:,i_eig)   !-Zim
      end do
!
!     compute the residuals, and their rms and sup norms:
!
!     real part 
!
      call dgemm('n','n',n,n_max,ldu,one,lvpre,n,up,lda,zero,rpre,n)  
      call dgemm('n','n',n,n_max,ldu,one,lvmre,n,um,lda,zero,rmre,n)
      call dgemm('n','n',n,n_max,ldu,one,bvpre,n,um,lda,zero,bpre,n)
      call dgemm('n','n',n,n_max,ldu,one,bvmre,n,up,lda,zero,bmre,n)
!
!     imaginary part 
!
      call dgemm('n','n',n,n_max,ldu,one,lvpim,n,up,lda,zero,rpim,n)
      call dgemm('n','n',n,n_max,ldu,one,lvmim,n,um,lda,zero,rmim,n)
      call dgemm('n','n',n,n_max,ldu,one,bvpim,n,um,lda,zero,bpim,n)
      call dgemm('n','n',n,n_max,ldu,one,bvmim,n,up,lda,zero,bmim,n)
      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
!       real part
!
        call daxpy(n,-eig(i_eig),bpre(:,i_eig),1,rpre(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bmre(:,i_eig),1,rmre(:,i_eig),1)
!
!       imaginary part
!
        call daxpy(n,-eig(i_eig),bpim(:,i_eig),1,rpim(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bmim(:,i_eig),1,rmim(:,i_eig),1)
!        
        r_norm(1,i_eig) = sqrt(dot_product(rpre(:,i_eig),rpre(:,i_eig)) + &
                               dot_product(rmre(:,i_eig),rmre(:,i_eig)) + &
                               dot_product(rpim(:,i_eig),rpim(:,i_eig)) + &
                               dot_product(rmim(:,i_eig),rmim(:,i_eig)))/sqrt2n
        r_norm(2,i_eig) = maxval(abs(rpre(:,i_eig) + rmre(:,i_eig) + & 
                                     rpim(:,i_eig) + rmim(:,i_eig)))
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
        done(i_eig+1)   = r_norm(1,i_eig+1).lt.tol_rms .and. &
                          r_norm(2,i_eig+1).lt.tol_max .and. &
                          it.gt.1
!                  
!     check coupled-eigenvalues convergece: if one eigval from the same couple
!     is false, both need to be put false (even if the other was .true.)
!
        if (done(i_eig) .neqv. done(i_eig+1)) then 
          done(i_eig)   = .false.
          done(i_eig+1) = .false.
        end if
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
!
!       real part 
!
        call lrprec(n,n_act,eig(ind),rpre(1,ind),rmre(1,ind),vpre(1,i_beg),vmre(1,i_beg))
!
!       imaginary part 
!
        call lrprec(n,n_act,eig(ind),rpim(1,ind),rmim(1,ind),vpim(1,i_beg),vmim(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
!        
        vpp(1:n,:)    = vpre
        vpp(n+1:n2,:) = vpim
        vmm(1:n,:)    = vmre
        vmm(n+1:n2,:) = vmim
        call ortho_vs_x_newcomplex(.false.,n2,ldu,n_act,vpp,vpp(1,i_beg),xx,xx)
        call ortho_vs_x_newcomplex(.false.,n2,ldu,n_act,vmm,vmm(1,i_beg),xx,xx)
        vpre = vpp(1:n,:)
        vpim = vpp(n+1:n2,:)
        vmre = vmm(1:n,:)
        vmim = vmm(n+1:n2,:)
!
!       real part 
!
!        call ortho_vs_x(.false.,n,ldu,n_act,vpre,vpre(1,i_beg),xx,xx)
!        call ortho_vs_x(.false.,n,ldu,n_act,vmre,vmre(1,i_beg),xx,xx)
!
!       imaginary part 
!
!        call ortho_vs_x(.false.,n,ldu,n_act,vpim,vpim(1,i_beg),xx,xx)
!        call ortho_vs_x(.false.,n,ldu,n_act,vmim,vmim(1,i_beg),xx,xx)
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
        vpre = zero
        vmre = zero
        vpim = zero
        vmim = zero
        vpp  = zero
        vmm  = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max, 2
!         real 
          vpre(:,i_eig)   = ( evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
          vmre(:,i_eig)   = ( evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
          vpre(:,i_eig+1) = ( evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
          vmre(:,i_eig+1) = ( evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!         imaginary 
          vpim(:,i_eig)   = ( evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
          vmim(:,i_eig)   = ( evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
          vpim(:,i_eig+1) = (-evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
          vmim(:,i_eig+1) = (-evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
        end do
!
        vpp(1:n,:)     = vpre
        vpp(n+1:n2,:)  = vpim
        vmm(1:n,:)     = vmre
        vmm(n+1:n2,:)  = vmim
        call ortho_cd_newcomplex(.false.,n2,n_max,vpp,xx,ok)
        call ortho_cd_newcomplex(.false.,n2,n_max,vmm,xx,ok)
        vpre = vpp(1:n,:)
        vpim = vpp(n+1:n2,:)
        vmre = vmm(1:n,:)
        vmim = vmm(n+1:n2,:)
!
!       real 
!        call ortho_cd(.false.,n,n_max,vpre,xx,ok)
!        call ortho_cd(.false.,n,n_max,vmre,xx,ok)
!       imaginary
!        call ortho_cd(.false.,n,n_max,vpim,xx,ok)
!        call ortho_cd(.false.,n,n_max,vmim,xx,ok)
!
        lvpre   = zero
        lvpim   = zero
        lvmre   = zero
        lvmim   = zero
        bvpre   = zero
        bvpim   = zero
        bvmre   = zero
        bvmim   = zero
        a_red   = zero
        s_red   = zero
        epmatre = zero
        epmatim = zero
        epmat   = zero
        emmatre = zero
        emmatim = zero
        emmat   = zero
        smatre  = zero
        smatim  = zero
        smat    = zero
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
    deallocate(work,tau,vpre,vmre,lvpre,lvmre,bvpre,bvmre,rpre,rmre,rrre, & 
               vpim,vmim,lvpim,lvmim,bvpim,bvmim,rpim,rmim,rrim,done,r_norm, & 
               a_red,a_copy,s_red,s_copy,e_red,epmatre,epmatim,emmatre,emmatim, &
               smatre,smatim,epmat,emmat,smat,up,um,eigpre,eigmre,eigpim,eigmim, & 
               vpp,vmm,bpre,bmre,bpim,bmim)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_newcomplex_driver
!
  subroutine caslr_complex_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                                  apbmul_complex,ambmul_complex,spdmul_complex, &
                                  smdmul_complex,lrprec_complex,eig,ceig,evec,abbamat,sddsmat,ok)
    use utils
    implicit none
    logical,                          intent(in)    :: verbose
    integer,                          intent(in)    :: n, n2, n_targ, n_max
    integer,                          intent(in)    :: max_iter, max_dav
    real(dp),                         intent(in)    :: tol
    real(dp),    dimension(n_max),    intent(inout) :: eig
    !real(dp),    allocatable                        :: rwork(:)
    !complex(dp), allocatable                        :: cwork(:)
    complex(dp), dimension(n_max),    intent(inout) :: ceig
    complex(dp), dimension(n2,n_max), intent(inout) :: evec
    complex(dp), dimension(n2,n2),    intent(inout) :: abbamat, sddsmat
    logical,                          intent(inout) :: ok
    external                                        :: apbmul_complex, ambmul_complex, spdmul_complex, smdmul_complex, &
                                                       lrprec_complex
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
    real(dp)              :: sqrtn, sqrt2n, tol_rms, tol_max
    complex(dp)           :: xx(1)
    real(dp)              :: lw_svd(1), lw(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    complex(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:), vpm(:,:), prodscr(:,:)
    complex(dp), allocatable :: rp(:,:), rm(:,:), rr(:,:)
    real(dp),    allocatable :: r_norm(:,:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp),    allocatable :: up(:,:), um(:,:)
    complex(dp), allocatable :: cup(:,:), cum(:,:), eigp(:,:), eigm(:,:), bp(:,:), bm(:,:)
    complex(dp), allocatable :: eigp_(:,:), eigm_(:,:)
    complex(dp), allocatable :: bppiu(:,:), bmmeno(:,:) !debug matrices b+, b- 2nx2n 
!
!   subspace matrix and eigenvalues.
!
    real(dp),    allocatable :: epmat(:,:), emmat(:,:), smat(:,:), e_red(:), eigen_a(:)
    complex(dp), allocatable :: a_red(:,:), a_copy(:,:), a_hlf(:,:), s_red(:,:), s_copy(:,:), cscr2(:,:)
    complex(dp), allocatable :: csppmat(:,:), csmmmat(:,:), cspmmat(:,:), csmpmat(:,:)
    complex(dp), allocatable :: eppiu(:,:), emmeno(:,:), esse(:,:), cscr(:,:)  !debug matrices 2nx2n
    complex(dp), allocatable :: cepmat(:,:), cemmat(:,:), csmat(:,:)!dc, s_red(:,:), s_copy(:,:), e_red(:)
    complex(dp), allocatable :: cepmmat(:,:), cempmat(:,:), supererre(:,:), minorerre(:,:)
!
!   to debug
!
    real(dp), allocatable    :: vpreal(:,:), vp_real(:,:), vp_imag(:,:), vm_real(:,:), vm_imag(:,:)
    integer                  :: lrwork, liwork
    integer, allocatable     :: iwork(:)
    complex(dp), allocatable :: evec_test(:,:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2, dznrm2
    external              :: dcopy, dnrm2, dgemm, dsyev, zcopy, dznrm2, zgemm, zheevd, zhegv
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    write(6,*) ' TRAD-COMPLEX OLSEN Diaglib '
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
    lda2    = 2 * lda
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), tau(n_max), ctau(n_max), cwork(lwork), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    allocate (vp(n,lda), vm(n,lda), lvp(n,lda), lvm(n,lda), bvp(n,lda), bvm(n,lda), &
              rp(n,n_max), rm(n,n_max), rr(n,n_max), supererre(n2,n_max), minorerre(n,n_max), &
              bppiu(2*n,lda), bmmeno(2*n,lda), vpm(n,lda2), prodscr(n,lda2), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (a_red(lda2,lda2), a_copy(lda2,lda2), a_hlf(lda2,lda), s_red(lda2,lda2), s_copy(lda2,lda2), &
              e_red(lda2), epmat(lda,lda), cepmat(lda,lda), emmat(lda,lda), cemmat(lda,lda), &
              cepmmat(lda,lda), cempmat(lda,lda), csppmat(lda,lda), csmmmat(lda,lda), &
              smat(lda,lda), csmat(lda,lda), eppiu(lda,lda), cscr(n2,n2), emmeno(lda,lda), esse(lda,lda), &
              evec_test(n2,n_max), cscr2(lda2,lda2), cspmmat(lda,lda), csmpmat(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), cup(lda,n_max), cum(lda,n_max), eigp(n,n_max), &
              eigm(n,n_max), bp(n,n_max), bm(n,n_max), eigp_(n,n_max), eigm_(n,n_max), stat = istat)
    call check_mem(istat)
!dc test
    allocate(vpreal(n,lda))
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn    = sqrt(real(n,dp))
    tol_rms  = tol
    tol_max  = 10.0_dp * tol
!
!   clean out various quantities
!
    t_diag  = zero
    t_ortho = zero
    t_mv    = zero
    t_tot   = zero
    vp      = czero
    vm      = czero
    vpm     = czero
    prodscr = czero
    bvp     = czero
    bvm     = czero
    lvp     = czero
    lvm     = czero
    cup     = czero
    cum     = czero
    !up      = zero 
    !um      = zero 
    a_red   = czero 
    s_copy  = czero
    cscr    = czero
    supererre = czero
    minorerre = czero
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
!    write(6,*) 'evec', evec(1:n,1)
!    write(6,*) 'evec2', evec(n+1:n2,1)
    do i_eig = 1, n_max
      vp(:,i_eig) = evec(1:n,i_eig) + evec(n+1:n2,i_eig)
      vm(:,i_eig) = evec(1:n,i_eig) - evec(n+1:n2,i_eig)
      !vp(:,i_eig) = evec(1:n,i_eig) + conjg(evec(n+1:n2,i_eig))
      !vm(:,i_eig) = evec(1:n,i_eig) - conjg(evec(n+1:n2,i_eig))
    end do
!    write(6,'(2a15)') 'vp','vm'
!    write(6,'(2f15.8)') real(vp), real(vm)
    !write(6,*) 'vp', real(vp)
    !write(6,*) 'vpim', aimag(vp)
    !write(6,*) 'vm', real(vm)
    !write(6,*) 'vmim', aimag(vm)
!    stop
!    print *, 'v +:'
!    vpreal = vp
!    print '(10f15.8)', transpose(vpreal(1:n_max,1:n_max)) 
    !print *, real(vm(:,1))
    !stop
!
!    allocate(vp_real(n,lda),vp_imag(n,lda),vm_real(n,lda),vm_imag(n,lda))
!    write(6,*) 'before ortho CD'
!    write(6,'(4a15)') 'Re vp','Im vp', 'Re vm', 'Im vm'
!    vp_real = zero
!    vm_real = zero
!    vp_imag = zero
!    vm_imag = zero
!    do j = 1, lda
!      do i = 1, n
!        vp_real(i,j) = real(vp(i,j))
!        vm_real(i,j) = real(vm(i,j))
!        vp_imag(i,j) = dimag(vp(i,j))
!        vm_imag(i,j) = dimag(vm(i,j))
!      end do
!    end do
!    do i = 1, n
!      write(6,'(4f15.8)') vp_real(i,1), vp_imag(i,1), vm_real(i,1), vm_imag(i,1)
!    end do
    call ortho_cd_complex(.false.,n,n_max,vp,xx,ok)
    call ortho_cd_complex(.false.,n,n_max,vm,xx,ok)
!    write(6,*) 'xx', xx
    !write(6,*) 'vp times vm:', dot_product(vp(:,1), vm(:,1))
!    write(6,*) 'after ortho CD'
!    write(6,'(4a15)') 'Re vp','Im vp', 'Re vm', 'Im vm'
!    do j = 1, lda
!      do i = 1, n
!        vp_real(i,j) = real(vp(i,j))
!        vm_real(i,j) = real(vm(i,j))
!        vp_imag(i,j) = dimag(vp(i,j))
!        vm_imag(i,j) = dimag(vm(i,j))
!      end do
!    end do
!    do i = 1, n
!      write(6,'(4f15.8)') vp_real(i,1), vp_imag(i,1), vm_real(i,1), vm_imag(i,1)
!    end do
!    deallocate(vp_real, vp_imag, vm_real, vm_imag)
!    stop
!
    call ortho_vs_x_complex(.false.,n,n_max,n_max,vp,vm,xx,xx,0)
    !write(6,*) 'vp scalare vm'
    !print *, dot_product(vp(:,1),vm(:,1)) - conjg(dot_product(vp(:,1),vm(:,1)))
    allocate(vp_real(n,lda),vp_imag(n,lda),vm_real(n,lda),vm_imag(n,lda))
    write(6,*) 'after ortho vs x'
    write(6,'(4a15)') 'Re vp','Im vp', 'Re vm', 'Im vm'
    vp_real = zero
    vm_real = zero
    vp_imag = zero
    vm_imag = zero
    do j = 1, lda
      do i = 1, n
        vp_real(i,j) = real(vp(i,j))
        vm_real(i,j) = real(vm(i,j))
        vp_imag(i,j) = dimag(vp(i,j))
        vm_imag(i,j) = dimag(vm(i,j))
      end do
    end do
    do i = 1, n
      write(6,'(4f15.8)') vp_real(i,1), vp_imag(i,1), vm_real(i,1), vm_imag(i,1)
    end do
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
    1030 format(t5,'Traditional complex driver iterations (tol=',d10.2,'):',/, &
            t5,'------------------------------------------------------------------',/, &
            t7,'  iter  root              eigenvalue','         rms         max ok',/, &
            t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    allocate(eigen_a(lda2))    !to debug, remove when it works 
    n_rst   = 0
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ldu = ldu + n_act
!
!     to verify whether the vectors are orthogonal
!
      !dccall zcopy(n*ldu,vp,1,vpm,1)
      !dccall zcopy(n*ldu,vm,1,vpm(1,ldu+1),1)
      !dcprodscr = czero
      !dccall zgemm('c','n',ldu,ldu,n,cone,vpm,n,vpm,n,czero,prodscr,n)   !If orthogonal, prodscr = identity 
      !dcprodscr = prodscr + conjg(prodscr)
      !dcwrite(6,*) 'prodscr' 
      !dcdo i = 1, ldu
      !dc  do j = 1, ldu
      !dc    write(6,*) 'i, j = ', i, j
      !dc    write(6,*) prodscr(i,j)
      !dc  end do
      !dcend do
!
!     perform this iteration's matrix-vector multiplication:
!
      call get_time(t1)
      call apbmul_complex(n,n_act,vp(1,i_beg),lvp(1,i_beg))
      call ambmul_complex(n,n_act,vm(1,i_beg),lvm(1,i_beg))
!      write(6,*) 'vp pre spd', vp
      call spdmul_complex(n,n_act,vp(1,i_beg),bvm(1,i_beg))
!      write(6,*) 'bvm post spd', bvm
!      write(6,*) 'vm pre smd', vm
      call smdmul_complex(n,n_act,vm(1,i_beg),bvp(1,i_beg))
!      write(6,*) 'bvp post smd', bvp
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     trying to debug the code: build the whole E+,E- matrices
!
!     b(+) = (b+ b+*)      b(-) = (b- -b-*)
!      do i = 1, ldu  
!        bppiu(1:n,i)      = vp(:,i)
!        bppiu(n+1:2*n,i)  = conjg(vp(:,i))
!        bmmeno(1:n,i)     = vm(:,i)
!        bmmeno(n+1:2*n,i) = -conjg(vm(:,i))
!      end do
!      eppiu  = czero
!      emmeno = czero 
!      esse   = czero
!      cscr   = czero 
!      !write(6,*) 'abbamat', abbamat !Printing mat AB,B*A* to debug
!      !write(6,*) 'sddsmat', sddsmat !Printing mat SD,-S*-D* to debug
!      call zgemm('n','n',n2,ldu,n2,cone,abbamat,n2,bppiu,n2,czero,cscr,n2)
!      call zgemm('c','n',ldu,ldu,n2,cone,bppiu,n2,cscr,n2,czero,eppiu,lda)
!      write(6,*) 'E+ da zgemm'
!      epmat = aimag(eppiu)
!      !write(6,'(10f20.10)') transpose(epmat(1:10,1:10))
!      write(6,'(f20.10)') maxval(abs(aimag(eppiu)))
!      call zgemm('n','n',n2,ldu,n2,cone,abbamat,n2,bmmeno,n2,czero,cscr,n2)
!      call zgemm('c','n',ldu,ldu,n2,cone,bmmeno,n2,cscr,n2,czero,emmeno,lda)
!      write(6,*) 'E- da zgemm'
!      emmat = aimag(emmeno)
!      write(6,'(f20.10)') maxval(abs(aimag(emmeno)))
!      call zgemm('n','n',n2,ldu,n2,cone,sddsmat,n2,bppiu,n2,czero,cscr,n2)
!      call zgemm('c','n',ldu,ldu,n2,cone,bmmeno,n2,cscr,n2,czero,esse,lda)
!      write(6,*) 'S da zgemm'
!      smat = aimag(esse)
!      !write(6,'(10f20.10)') transpose(smat(1:10,1:10))
!      write(6,'(f20.10)') maxval(abs(aimag(esse)))
!
!     update the reduced matrix 
!
!     build the blocks E++,E--,E+-,E-+
!
      call zgemm('c','n',ldu,ldu,n,cone,vp,n,lvp,n,czero,cepmat,lda)   !E++
      cepmat = cepmat + conjg(cepmat)
      !write(6,*) 'epmat', epmat
      call zgemm('c','n',ldu,ldu,n,cone,vp,n,lvm,n,czero,cepmmat,lda)  !E+-
      cepmmat = cepmmat - conjg(cepmmat)
      !write(6,*) 'cepmmat', cepmmat
      call zgemm('c','n',ldu,ldu,n,cone,vm,n,lvp,n,czero,cempmat,lda)  !E-+
      cempmat = cempmat - conjg(cempmat)
!      cempmat = transpose(conjg(cepmmat)) 
      !write(6,*) 'cempmat', cempmat
      call zgemm('c','n',ldu,ldu,n,cone,vm,n,lvm,n,czero,cemmat,lda)   !E--
      cemmat = cemmat + conjg(cemmat)
      !write(6,*) 'emmat', emmat
!
!    build the blocks S++,S,S-- 
!
      call zgemm('c','n',ldu,ldu,n,cone,vp,n,bvm,n,czero,csppmat,lda)    !S++
      csppmat  = csppmat - conjg(csppmat)
      write(6,*) "S++"
      do i = 1, ldu
        do j = 1, ldu
          write(6,*) csppmat(j,i)
        end do
      end do
!      write(6,*) 'csppmat'
!      write(6,'(2f10.8)') real(csppmat), aimag(csppmat)
      !call zgemm('c','n',ldu,ldu,n,cone,vm,n,bvm,n,czero,csmat,lda)    !S-+=S
      !smat  = csmat + conjg(csmat)
      call zgemm('c','n',ldu,ldu,n,cone,vp,n,bvp,n,czero,cspmmat,lda)    !S+-
      cspmmat  = cspmmat + conjg(cspmmat)
      write(6,*) "S+-"
      do i = 1, ldu
        do j = 1, ldu
          write(6,*) cspmmat(j,i)
        end do
      end do
      call zgemm('c','n',ldu,ldu,n,cone,vm,n,bvm,n,czero,csmpmat,lda)    !S-+
      csmpmat  = csmpmat + conjg(csmpmat)
      !write(6,*) 'smat', smat
      call zgemm('c','n',ldu,ldu,n,cone,vm,n,bvp,n,czero,csmmmat,lda)    !S--
      csmmmat  = csmmmat - conjg(csmmmat)
      !write(6,*) 'csmmmat', csmmmat
!
      !smat  = real(csmat)
      !epmat = real(cepmat) 
      !emmat = real(cemmat) 
      !write(6,*) 'cE+, f.pes portami a cagliari con te'
      !write(6,*) cepmat
      !write(6,*) 'cE-, f.pes portami a cagliari con te'
      !write(6,*) cemmat
      !write(6,*) 'cS, f.pes portami a cagliari con te'
      !write(6,*) csmat
      !smat  = csmat + conjg(csmat)
      !write(6,*) 'S=cS + conjg(cS)'
      !write(6,'(10f20.10)') transpose(smat(1:10,1:10))
      !epmat = cepmat + conjg(cepmat)
      !write(6,*) 'E+=cE+ + conjg(cE+)'
      !write(6,'(10f20.10)') transpose(epmat(1:10,1:10))
      !emmat = cemmat + conjg(cemmat)
      !write(6,*) 'E-=cE- + conjg(cE-)'
      !write(6,'(10f20.10)') transpose(emmat(1:10,1:10))
      !smat  = two*real(csmat)
      !epmat = two*real(cepmat)
      !emmat = two*real(cemmat)
      a_red = czero 
      s_red = czero
      !a_red(1:ldu,1:ldu)             = epmat(1:ldu,1:ldu)
      a_red(1:ldu,1:ldu)             = cepmat(1:ldu,1:ldu)
      a_red(1:ldu,ldu+1:2*ldu)       = cepmmat(1:ldu,1:ldu)
!      a_red(ldu+1:2*ldu,1:ldu)       = transpose(conjg(cepmmat(1:ldu,1:ldu)))  !quando funzionerà tieni solo questo
      a_red(ldu+1:2*ldu,1:ldu)       = cempmat(1:ldu,1:ldu)
      !a_red(ldu+1:2*ldu,ldu+1:2*ldu) = emmat(1:ldu,1:ldu)
      a_red(ldu+1:2*ldu,ldu+1:2*ldu) = cemmat(1:ldu,1:ldu)
      s_red(1:ldu,1:ldu)             = csppmat(1:ldu,1:ldu)
      !s_red(ldu+1:2*ldu,1:ldu)       = smat(1:ldu,1:ldu)
      !s_red(1:ldu,ldu+1:2*ldu)       = transpose(smat(1:ldu,1:ldu))
      s_red(ldu+1:2*ldu,1:ldu)       = csmpmat(1:ldu,1:ldu)
      s_red(1:ldu,ldu+1:2*ldu)       = cspmmat(1:ldu,1:ldu)
      s_red(ldu+1:2*ldu,ldu+1:2*ldu) = csmmmat(1:ldu,1:ldu)
!      write(6,*) 's_red'
!      call prtmat_complex(2*ldu,2*ldu,s_red,3)
!      stop
!
!     to verify that a_red is positive definite
!
!      e_red = zero
!      a_copy = a_red
!      deallocate(cwork)
!      allocate (cwork(1), rwork(1), iwork(1))
!      call zheevd('v','l',2*ldu,a_copy,2*lda,eigen_a,cwork,-1,rwork,-1,iwork,-1,info)
!      write(6,*) 'info1', info
!      lwork = int(cwork(1))
!      lrwork = int(rwork(1))
!      liwork = int(iwork(1))
!      deallocate(cwork, rwork, iwork)
!      allocate (cwork(lwork), rwork(lrwork), iwork(liwork))
!      call zheevd('v','l',2*ldu,a_copy,2*lda,eigen_a,cwork,lwork,rwork,lrwork,iwork,liwork,info)
!      deallocate(rwork,iwork)
!      write(6,*) 'info2', info
!      j = 0
!      do i = 1, 2*ldu
!        if (eigen_a(i) .le. 0.0_dp) then
!          write(6,*) 'WARNING: NOT POSITIVE DEFINITE :(, eigenvalues:', eigen_a
!        else
!          j = j + 1
!        end if 
!      end do
!      if (j .eq. 2*ldu) write(6,*) 'a_red is positive definite, we are happy :)', eigen_a
!
!    to verify if a_red is hermitian
!
!      write(6,*) 'a_red'
!      write(6,*) 'if norm of a_red is zero, we are happy:)', dznrm2(2*ldu,a_red - transpose(conjg(a_red)),1)
!      write(6,*) 'a_red'
!      call prtmat_complex(2*ldu,2*ldu,a_red,3)
!
!    to verify if s_red is hermitian
!
!      write(6,*) 's_red'
!      call prtmat_complex(2*ldu,2*ldu,s_red,3)
!      write(6,*) 'if norm of s_red is zero, we are happy:)', dznrm2(2*ldu,s_red - transpose(conjg(s_red)),1)
!      stop
!
      call get_time(t1)
!
!     default algorithm: solve the 2n-dimensional inverse problem.
!
      deallocate(cwork)
      allocate(cwork(1), rwork(1))
      call zhegv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,cwork,-1,rwork,info)
      if (info .ne. 0) then 
        write(6,*) 'WARNING caslr_complex dryrun: info zhegv', info
        stop
      end if
      lwork = int(cwork(1))
      deallocate(cwork,rwork)
      allocate(cwork(10*lwork),rwork(300*(2*ldu)-2))
      call zhegv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,cwork,lwork,rwork,info)
      deallocate(rwork)
      if (info .ne. 0) then 
        write(6,*) 'WARNING caslr_complex: info zhegv', info
        stop
      end if
!
!    real diagonalization to solve the 2n-dimensional inverse problem (traditional algorithm)
!
!      deallocate(work)
!      call dsygv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,lw,-1,info)
!      if (info .ne. 0) then 
!        write(6,*) 'WARNING caslr_complex: info dsygv', info
!        stop
!      end if
!      lwork = int(lw(1))
!      allocate(work(lwork))
!      call dsygv(1,'v','l',2*ldu,s_red,2*lda,a_red,2*lda,e_red,work,lwork,info)
!      if (info .ne. 0) then 
!        write(6,*) 'WARNING caslr_complex: info dsygv', info
!        stop
!      end if
!
!     new diagonalization of s_red avoiding the lapack
!
!     first step: diagonalize a_red
!
!dc      e_red = zero
!dc      a_copy = a_red
!dc      deallocate(cwork)
!dc      allocate (cwork(1), rwork(1), iwork(1))
!dc      call zheevd('v','l',2*ldu,a_copy,2*lda,eigen_a,cwork,-1,rwork,-1,iwork,-1,info)
!dc      write(6,*) 'info1', info
!dc      lwork = int(cwork(1))
!dc      lrwork = int(rwork(1))
!dc      liwork = int(iwork(1))
!dc      deallocate(cwork, rwork, iwork)
!dc      allocate (cwork(lwork), rwork(lrwork), iwork(liwork))
!dc      call zheevd('v','l',2*ldu,a_copy,2*lda,eigen_a,cwork,lwork,rwork,lrwork,iwork,liwork,info)
!dc      deallocate(rwork,iwork)
!dc      write(6,*) 'info2', info
!dc!
!dc!     second step: compute eigenval**-1/2
!dc!
!dc      do i = 1, 2*ldu
!dc        if (abs(eigen_a(i)) .le. 1e-7) then 
!dc          write(6,*) 'eigenval number', i, 'is zero'
!dc          eigen_a(i) = zero
!dc        else 
!dc          eigen_a(i) = one/sqrt(eigen_a(i))
!dc        end if
!dc      end do
!dc      write(6,*) ' eigenvalues^-1/2'
!dc      write(6,'(10f15.8)') eigen_a(1:10)
!dc!
!dc!     third step: compute a_red**-1/2 = U*eigenval**-1/2*U**H
!dc!
!dc      cscr2 = czero
!dc      do i = 1, 2*ldu
!dc        cscr2(:,i) = a_copy(:,i)*dcmplx(eigen_a(i)) 
!dc      end do
!dc      call zgemm('n','c',2*ldu,2*ldu,2*ldu,cone,cscr2,lda2,a_copy,lda2,czero,a_hlf,lda2)
!dc!
!dc!     fourth step: compute s_red' = a_hlf*s_red_*a_hlf
!dc!    
!dc      call zgemm('n','n',2*ldu,2*ldu,2*ldu,cone,s_red,lda2,a_hlf,lda2,czero,a_copy,lda2)
!dc      call zgemm('n','n',2*ldu,2*ldu,2*ldu,cone,a_hlf,lda2,a_copy,lda2,czero,s_red,lda2)
!dc!
!dc!     fifth step: diagonalize s_red' --> s_red'x' = lambda x'
!dc!
!dc      s_copy = s_red
!dc      deallocate(cwork)
!dc      allocate (cwork(1), rwork(1), iwork(1))
!dc      call zheevd('v','l',2*ldu,s_copy,2*lda,e_red,cwork,-1,rwork,-1,iwork,-1,info)
!dc      write(6,*) 'info1', info
!dc      lwork = int(cwork(1))
!dc      lrwork = int(rwork(1))
!dc      liwork = int(iwork(1))
!dc      deallocate(cwork, rwork, iwork)
!dc      allocate (cwork(lwork), rwork(lrwork), iwork(liwork))
!dc      call zheevd('v','l',2*ldu,s_copy,2*lda,e_red,cwork,lwork,rwork,lrwork,iwork,liwork,info)
!dc      deallocate(rwork,iwork)
!dc      write(6,*) 'info2', info
!dc!
!dc!     last step: obtain the eigenvector x = a_hlf*x'
!dc!
!dc      call zgemm('n','n',2*ldu,2*ldu,2*ldu,cone,a_hlf,lda2,s_copy,lda2,czero,s_red,lda2)
!      write(6,*) 'x final', s_red
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
!      write(6,*) 'e_red', e_red
      do i_eig = 1, n_max
        eig(i_eig)       = one/e_red(2*ldu - i_eig + 1)
        !up(1:ldu,i_eig)  = s_red(1:ldu,2*ldu - i_eig + 1)
        !um(1:ldu,i_eig)  = s_red(ldu+1:2*ldu,2*ldu - i_eig + 1)
        cup(1:ldu,i_eig)  = s_red(1:ldu,2*ldu - i_eig + 1)
        cum(1:ldu,i_eig)  = s_red(ldu+1:2*ldu,2*ldu - i_eig + 1)
        !cup(1:ldu,i_eig) = up(1:ldu,i_eig)
        !cum(1:ldu,i_eig) = um(1:ldu,i_eig)
!fl
!       write(6,*) 'eig: ', eig(i_eig)
!       write(6,*) 'up ', i_eig
!       write(6,'(10f12.6)') up(1:ldu,i_eig)
!       write(6,*) 'um ', i_eig
!       write(5,'(10f12.6)') um(1:ldu,i_eig)
      end do
!
      !write(6,*) 'up', up
      !write(6,*) 'um', um
      !write(6,*) 'cup', cup
      !write(6,*) 'cum', cum
      call get_time(t2)
      t_diag = t_diag + t2 - t1
      call zgemm('n','n',n,n_max,ldu,cone,vp,n,cup,lda,czero,eigp,n)
      call zgemm('n','n',n,n_max,ldu,cone,vm,n,cum,lda,czero,eigm,n)
      call zgemm('n','n',n,n_max,ldu,cone,conjg(vp),n,cup,lda,czero,eigp_,n)
      call zgemm('n','n',n,n_max,ldu,cone,conjg(vm),n,cum,lda,czero,eigm_,n)
!      if (it.eq.3) print *, vp(:,2)
!      if (it.eq.3) stop
!
      !write(6,*) 'vp', vp
      !write(6,*) 'vp*', conjg(vp)
!      write(6,*) 'eigp', eigp
!      write(6,*) 'eigp_', eigp_
!      write(6,*) 'eigm', eigm
!      write(6,*) 'eigm_', eigm_
      do i_eig = 1, n_max
        evec(1:n,i_eig)    = eigp(:,i_eig) + eigm(:,i_eig)
!        evec(n+1:n2,i_eig) = eigp(:,i_eig) - eigm(:,i_eig)
        evec(n+1:n2,i_eig) = eigp_(:,i_eig) - eigm_(:,i_eig)
        !evec(n+1:n2,i_eig) = conjg(evec(n+1:n2,i_eig))
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call zgemm('n','n',n,n_max,ldu,cone,lvp,n,cup,lda,czero,rp,n)
      call zgemm('n','n',n,n_max,ldu,cone,lvm,n,cum,lda,czero,rm,n)
      call zgemm('n','n',n,n_max,ldu,cone,bvp,n,cum,lda,czero,bp,n)
      call zgemm('n','n',n,n_max,ldu,cone,bvm,n,cup,lda,czero,bm,n)
!      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        ceig(i_eig) = dcmplx(eig(i_eig))
        call zaxpy(n,-ceig(i_eig),bp(:,i_eig),1,rp(:,i_eig),1)
        call zaxpy(n,-ceig(i_eig),bm(:,i_eig),1,rm(:,i_eig),1)
!        r_norm(1,i_eig) = dznrm2(n,rp(:,i_eig),1)/sqrtn + dznrm2(n,rm(:,i_eig),1)/sqrtn
        r_norm(1,i_eig) = dznrm2(n,rp(:,i_eig) + rm(:,i_eig),1)/sqrtn
!        r_norm(2,i_eig) = maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig)))
        r_norm(2,i_eig) = maxval(abs(rp(:,i_eig) + rm(:,i_eig)))
      end do
!
!     to debug: build the residuals (THE FUCKING WHOLE RESIDUALS)
!     !R = (abbamat - ceig*sddsmat)*evec
!      cscr = abbamat
!      do i_eig = 1, n_targ
!        if (done(i_eig)) cycle
!        call zaxpy(n2**2,-ceig(i_eig),sddsmat,1,cscr,1)
!        call zgemm('n','n',n2,1,n2,cone,cscr,n2,evec(:,i_eig),n2,czero,supererre(:,i_eig),n2)
!        cscr = abbamat
!        write(6,*) 'norm(rp + rm)', dznrm2(n,rp(:,i_eig) + rm(:,i_eig),1)/sqrtn
!        write(6,*) 'norm(supererre)', dznrm2(n2,supererre(:,i_eig),1)/sqrt(dble(n2))
!        write(6,*) 'max(|rp+rm|)', maxval(abs(rp(:,i_eig) + rm(:,i_eig)))
!        write(6,*) 'max(|supererre|)', maxval(abs(supererre(:,i_eig)))
!      end do
      !evec = evec/sqrt(ctwo)
!      write(6,*) 'supererre completo', supererre 
!      
!     R = R+ + R-
!    
!      minorerre = rp + rm
!      write(6,*) 'R = r+ + r-', supererre 
!      write(6,'(5a15)') 'Re(super)', 'Re(minor)', 'Im(super)', 'Im(minor)', '1 eig'
!      do i = 1, n
!        write(6,'(4f15.8)') real(supererre(i,1)), real(minorerre(i,1)), aimag(supererre(i,1)), aimag(minorerre(i,1))
!      end do
!      write(6,'(5a15)') 'Re(super)', 'Re(minor)', 'Im(super)', 'Im(minor)', '5 eig'
!      do i = 1, n
!        write(6,'(4f15.8)') real(supererre(i,5)), real(minorerre(i,5)), aimag(supererre(i,5)), aimag(minorerre(i,5))
!      end do
!      write(6,'(5a15)') 'Re(super)', 'Re(minor)', 'Im(super)', 'Im(minor)', '10 eig'
!      do i = 1, n
!        write(6,'(4f15.8)') real(supererre(i,10)), real(minorerre(i,10)), aimag(supererre(i,10)), aimag(minorerre(i,10))
!      end do
!      if (r_norm(1,1) .lt. 1.e-10_dp) return
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
        call lrprec_complex(n,n_act,ceig(ind),rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        !call ortho_vs_x_complex(.false.,n,ldu,n_act,vp,vp(1,i_beg),xx,xx,4)
        !call ortho_vs_x_complex(.false.,n,ldu,n_act,vm,vm(1,i_beg),xx,xx,4)
        call zcopy(n*ldu,vp,1,vpm,1)
        call zcopy(n*ldu,vm,1,vpm(1,ldu+1),1)
!        do i = 1, ldu
!          write(6,'(i2,a)') i, 'th vector:'
!          write(6,'(4a15)') 'Re(VP)','Re(VPM)VP','Re(VM)','Re(VPM)'
!          do j = 1, 5
!            write(6,'(4f15.8)') real(vp(j,i)), real(vpm(j,i)), real(vm(j,i)), real(vpm(j,ldu+i))
!          end do
!          write(6,*)
!        end do
!        write(6,*) 'vpm pre', vpm
        call ortho_vs_x_complex(.false.,n,2*ldu,n_act,vpm,vp(1,i_beg),xx,xx,0)
!        write(6,*) 'vpm post ', vpm
        call zcopy(n*n_act,vp(1,i_beg),1,vpm(1,(2*ldu)+1),1)
        call ortho_vs_x_complex(.false.,n,2*ldu+n_act,n_act,vpm,vm(1,i_beg),xx,xx,0)
        !do i_eig = 1, n_act
         ! print *, 'hello'
         !print *, dot_product(vp(:,i_beg+i_eig-1),vm(:,i_beg+i_eig-1))
         !print *, dot_product(vm(:,i_beg+i_eig-1),vm(:,i_beg+i_eig-1))
         !print *, dot_product(vp(:,i_eig),vp(:,i_beg+i_eig-1)) + conjg(dot_product(vp(:,i_eig),vp(:,i_beg+i_eig-1)))
         ! print *, dot_product(vm(:,i_eig),vm(:,i_beg+i_eig-1)) + conjg(dot_product(vm(:,i_eig),vm(:,i_beg+i_eig-1)))
         ! print *, dot_product(vp(:,i_eig),vm(:,i_beg+i_eig-1))-conjg(dot_product(vp(:,i_eig),vm(:,i_beg+i_eig-1)))
        !end do 
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting Davidson.'
        restart = .true.
!
!       initialize indices back to their starting values 
!
        ldu   = 0
        i_beg = 1
        m_dim = 1
        n_rst = 0
!
        n_act = n_max 
        vp    = czero
        vm    = czero
!        vpm   = czero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = evec(1:n,i_eig) + conjg(evec(n+1:n2,i_eig))
          vm(:,i_eig) = evec(1:n,i_eig) - conjg(evec(n+1:n2,i_eig))
        end do
!
        call ortho_cd_complex(.false.,n,n_max,vp,xx,ok)
        call ortho_cd_complex(.false.,n,n_max,vm,xx,ok)
        call ortho_vs_x_complex(.false.,n,n_max,n_max,vp,vm,xx,xx,0)
!
        lvp     = czero
        lvm     = czero
        bvp     = czero
        bvm     = czero
        a_red   = czero
        s_copy  = czero
        s_red   = czero
        cscr2   = czero
        a_hlf   = czero
        epmat   = zero
        emmat   = zero
        smat    = zero
        !cscr    = czero
        csppmat = czero
        csmmmat = czero
        cepmmat = czero
        cempmat = czero
        cepmat  = czero
        cemmat  = czero
        csmat   = czero
!
      end if
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!
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
    deallocate(work,cwork,ctau,tau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm,a_red,a_copy,s_red,s_copy,e_red, &
           epmat,emmat,smat,up,um,eigp,eigm,bp,bm,cepmat,cemmat,cepmmat,cempmat,csppmat,csmmmat,csmat,cup,cum,vpm,csmpmat,cspmmat)
!
    1050 format(t5,'----------------------------------------',/,&
        t7,'# target vectors:    ',i4,/,&
        t7,'# new vectors added: ',i4,/,&
        t7,'# converged vectors: ',i4,/,&
        t5,'----------------------------------------')
    return
  end subroutine caslr_complex_driver
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
!      
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
subroutine caslr_newcomplex_eff_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                                       apbmul,ambmul,spdmul,smdmul,lrprec,eig,evecre,evecim,ok)
!                               
    use utils
    implicit none
    logical, intent(in)                          :: verbose
    integer,                       intent(in)    :: n, n2, n_targ, n_max
    integer,                       intent(in)    :: max_iter, max_dav
    real(dp),                      intent(in)    :: tol
    real(dp), dimension(n_max),    intent(inout) :: eig
    real(dp), dimension(n2,n_max), intent(inout) :: evecre, evecim
    logical,                       intent(inout) :: ok
    external                                     :: apbmul, ambmul, spdmul, smdmul, &
                                                    lrprec
!
!   local variables:
!   ================
!
    integer, parameter    :: min_dav = 10
    integer               :: istat
    integer               :: i, j, k
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
    real(dp)              :: sqrtn, sqrt2n, tol_rms, tol_max
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
!   real
    real(dp), allocatable :: vpre(:,:), vmre(:,:), lvpre(:,:), lvmre(:,:), bvpre(:,:), bvmre(:,:)
    real(dp), allocatable :: rpre(:,:), rmre(:,:), rrre(:,:), r_norm(:,:)
!   imaginary
    real(dp), allocatable :: vpim(:,:), vmim(:,:), lvpim(:,:), lvmim(:,:), bvpim(:,:), bvmim(:,:)
    real(dp), allocatable :: rpim(:,:), rmim(:,:), rrim(:,:)
!
!   debug
!
    real(dp), allocatable :: vpp(:,:), vmm(:,:) 
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
!   real
    real(dp), allocatable :: up(:,:), um(:,:), eigpre(:,:), eigmre(:,:), bpre(:,:), bmre(:,:)
    real(dp), allocatable :: eigpre_(:,:), eigmre_(:,:)
!   imaginary
    real(dp), allocatable :: eigpim(:,:), eigmim(:,:), bpim(:,:), bmim(:,:)
    real(dp), allocatable :: eigpim_(:,:), eigmim_(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: smat(:,:), smatre(:,:), smatim(:,:)
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
    allocate (vpre(n,lda), vpim(n,lda), vmre(n,lda), vmim(n,lda), lvpre(n,lda), & 
              lvpim(n,lda), lvmre(n,lda), lvmim(n,lda), bvpre(n,lda), bvpim(n,lda), & 
              bvmre(n,lda), bvmim(n,lda), rpre(n,n_max), rpim(n,n_max), rmre(n,n_max), & 
              rmim(n,n_max), rrre(n,n_max), rrim(n,n_max), stat = istat)
    call check_mem(istat)
!
!   allocate memory for convergence check
!
    allocate (done(n_max), r_norm(2,n_max), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    allocate (s_red(lda,lda), s_copy(lda,lda), e_red(lda2), smat(lda,lda), & 
              smatre(lda,lda), smatim(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigpre(n,n_max), eigpim(n,n_max), &
              eigmre(n,n_max), eigmim(n,n_max), bpre(n,n_max), bpim(n,n_max), &
              bmre(n,n_max), bmim(n,n_max), eigpre_(n,n_max), eigmre_(n,n_max), & 
              eigpim_(n,n_max), eigmim_(n,n_max),stat = istat)
    call check_mem(istat)
!
!   debug
!
    allocate (vpp(n2,lda), vmm(n2,lda), stat = istat)
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    sqrt2n  = sqrt(real(n2,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
!
!   print something to understand we here in SMOGD
!
    write(6,*) 
    write(6,*) 'COMPLEX SMO-GD WITH REAL ALGEBRA IS ENTERED'
    write(6,*) 
!
!   clean out various quantities
!
    t_diag  = zero
    t_ortho = zero
    t_mv    = zero
    t_tot   = zero
    vpre    = zero
    vpim    = zero
    vmre    = zero
    vmim    = zero
    bvpre   = zero
    bvpim   = zero
    bvmre   = zero
    bvmim   = zero
    lvpre   = zero
    lvpim   = zero
    lvmre   = zero
    lvmim   = zero
    smat    = zero
    ok      = .false.
    done    = .false.
!
!   debug
!
    vpp = zero
    vmm = zero
!    
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max, 2
!     real 
      vpre(:,i_eig)   = (evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
      vmre(:,i_eig)   = (evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
      vpre(:,i_eig+1) = (-evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vmre(:,i_eig+1) = (-evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!     imaginary 
      vpim(:,i_eig)   = (evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig)   = (evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
      vpim(:,i_eig+1) = (evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig+1) = (evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
    end do
!    do i_eig = 1, n_max
!      write(6,*) 'vp 1 = r+ r+ i+ -i+', i_eig
!      do i = 1, n
!        write(6,'(4f8.4)') vpre(i,i_eig), vpre(i,i_eig), vpim(i,i_eig), -vpim(i,i_eig) 
!      end do
!      write(6,*) 'vm 2 = -i+ i+ r+ r+', i_eig
!      do i = 1, n
!        write(6,'(4f8.4)') vmre(i,i_eig), -vmre(i,i_eig), vmim(i,i_eig), vmim(i,i_eig) 
!      end do
!    end do 
!    stop
!
!   orthogonalize the expansion space to the metric.
!
    call apbmul(0,n,n_max,vpre,vpim,lvpre)
    call apbmul(1,n,n_max,vpre,vpim,lvpim)
    !
    call b_ortho(n,n_max,vpre,lvpre)
    call b_ortho(n,n_max,vpim,lvpim)
    !
    call ambmul(0,n,n_max,vmre,vmim,lvmre)
    call ambmul(1,n,n_max,vmre,vmim,lvmim)
    !
    call b_ortho(n,n_max,vmre,lvmre)
    call b_ortho(n,n_max,vmim,lvmim)
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
      call spdmul(0,n,n_act,vpre(1,i_beg),vpim(1,i_beg),bvmre(1,i_beg))
      call spdmul(1,n,n_act,vpre(1,i_beg),vpim(1,i_beg),bvmre(1,i_beg))
      call spdmul(0,n,n_act,vmre(1,i_beg),vmim(1,i_beg),bvpre(1,i_beg))
      call spdmul(1,n,n_act,vmre(1,i_beg),vmim(1,i_beg),bvpre(1,i_beg))
!
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvmre,n,zero,smatre,lda)
      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvmim,n,zero,smatim,lda)
      smat = smatre + smatim
      !smat = two*(smatre + smatim)
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
! TODO 
      call dgemm('n','n',n,n_max,ldu,one,vpre,n,up,lda,zero,eigpre,n)
      call dgemm('n','n',n,n_max,ldu,one,vmre,n,um,lda,zero,eigmre,n)
!
!     assemble the current approximation to the eigenvectors
!
! TODO 
      do i_eig = 1, n_max
        evecre(1:n,i_eig)    = eigpre(:,i_eig) + eigmre(:,i_eig)
        evecim(n+1:n2,i_eig) = eigpre(:,i_eig) - eigmre(:,i_eig)
      end do
!
!     compute the residuals, and their rms and sup norms:
!    
!     real part
!
      call dgemm('n','n',n,n_max,ldu,one,bvpre,n,um,lda,zero,rpre,n)
      call dgemm('n','n',n,n_max,ldu,one,bvmre,n,up,lda,zero,rmre,n)
      call dgemm('n','n',n,n_max,ldu,one,lvpre,n,up,lda,zero,bpre,n)
      call dgemm('n','n',n,n_max,ldu,one,lvmre,n,um,lda,zero,bmre,n)
!    
!     imaginary part
!
      call dgemm('n','n',n,n_max,ldu,one,bvpim,n,um,lda,zero,rpim,n)
      call dgemm('n','n',n,n_max,ldu,one,bvmim,n,up,lda,zero,rmim,n)
      call dgemm('n','n',n,n_max,ldu,one,lvpim,n,up,lda,zero,bpim,n)
      call dgemm('n','n',n,n_max,ldu,one,lvmim,n,um,lda,zero,bmim,n)
!      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
!       real part 
!
        call daxpy(n,-eig(i_eig),bpre(:,i_eig),1,rpre(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bmre(:,i_eig),1,rmre(:,i_eig),1)
!
!       imaginary part 
!
        call daxpy(n,-eig(i_eig),bpim(:,i_eig),1,rpim(:,i_eig),1)
        call daxpy(n,-eig(i_eig),bmim(:,i_eig),1,rmim(:,i_eig),1)
!        
!old        r_norm(1,i_eig) = (dnrm2(n,rp(:,i_eig),1) + dnrm2(n,rm(:,i_eig),1))/(eig(i_eig)*sqrt(two)*sqrtn)
!old        r_norm(2,i_eig) = (maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig))))/(sqrt(two)*eig(i_eig))
        r_norm(1,i_eig) = sqrt(dot_product(rpre(:,i_eig),rpre(:,i_eig)) + &
                               dot_product(rmre(:,i_eig),rmre(:,i_eig)) + &
                               dot_product(rpim(:,i_eig),rpim(:,i_eig)) + &
                               dot_product(rmim(:,i_eig),rmim(:,i_eig)))/sqrt2n
        r_norm(2,i_eig) = maxval(abs(rpre(:,i_eig) + rmre(:,i_eig) + & 
                                     rpim(:,i_eig) + rmim(:,i_eig)))
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
        done(i_eig+1)   = r_norm(1,i_eig+1).lt.tol_rms .and. &
                          r_norm(2,i_eig+1).lt.tol_max .and. &
                          it.gt.1
!                  
!     check coupled-eigenvalues convergece: if one eigval from the same couple
!     is false, both need to be put false (even if the other was .true.)
!
        if (done(i_eig) .neqv. done(i_eig+1)) then 
          done(i_eig)   = .false.
          done(i_eig+1) = .false.
        end if
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
!
!       real part 
!
        call lrprec(n,n_act,eig(ind),rpre(1,ind),rmre(1,ind),vpre(1,i_beg),vmre(1,i_beg))
!
!       imaginary part 
!
        call lrprec(n,n_act,eig(ind),rpim(1,ind),rmim(1,ind),vpim(1,i_beg),vmim(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        !
        call b_ortho_vs_x(n,ldu,n_act,vpre,lvpre,vpre(1,i_beg))
        call b_ortho_vs_x(n,ldu,n_act,vpim,lvpim,vpim(1,i_beg))
        !
        call apbmul(0,n,n_max,vpre,vpim,lvpre)
        call apbmul(1,n,n_max,vpre,vpim,lvpim)
        !
        call b_ortho(n,n_act,vpre(1,i_beg),lvpre(1,i_beg))
        call b_ortho(n,n_act,vpim(1,i_beg),lvpim(1,i_beg))
        !
        call b_ortho_vs_x(n,ldu,n_act,vmre,lvmre,vmre(1,i_beg))
        call b_ortho_vs_x(n,ldu,n_act,vmim,lvmim,vmim(1,i_beg))
        !
        call ambmul(0,n,n_max,vmre,vmim,lvmre)
        call ambmul(1,n,n_max,vmre,vmim,lvmim)
        !
        call b_ortho(n,n_act,vmre(1,i_beg),lvmre(1,i_beg))
        call b_ortho(n,n_act,vmim(1,i_beg),lvmim(1,i_beg))
        !
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
! TODO 
        n_act = n_max 
        vpre = zero
        vmre = zero
        vpim = zero
        vmim = zero
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max, 2
!         real 
          vpre(:,i_eig)   = (evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
          vmre(:,i_eig)   = (evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
          vpre(:,i_eig+1) = (-evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
          vmre(:,i_eig+1) = (-evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!         imaginary 
          vpim(:,i_eig)   = (evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
          vmim(:,i_eig)   = (evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
          vpim(:,i_eig+1) = (evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
          vmim(:,i_eig+1) = (evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
        end do
!
        lvpre   = zero
        lvmre   = zero
        lvpim   = zero
        lvmim   = zero
!
!       orthogonalize the expansion space to the metric.
!    
        call apbmul(0,n,n_max,vpre,vpim,lvpre)
        call apbmul(1,n,n_max,vpre,vpim,lvpim)
        !
        call b_ortho(n,n_max,vpre,lvpre)
        call b_ortho(n,n_max,vpim,lvpim)
        !
        call ambmul(0,n,n_max,vmre,vmim,lvmre)
        call ambmul(1,n,n_max,vmre,vmim,lvmim)
        !
        call b_ortho(n,n_max,vmre,lvmre)
        call b_ortho(n,n_max,vmim,lvmim)
!
        bvpre   = zero
        bvmre   = zero
        bvpim   = zero
        bvmim   = zero
        s_red   = zero
        smatre  = zero
        smatim  = zero
        smat    = zero
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
    deallocate(work,tau,vpre,vpim,vmre,vmim,lvpre,lvpim,lvmre,lvmim,bvpre,bvpim,bvmre,bvmim, &
               rpre,rpim,rmre,rmim,rrre,rrim,done,r_norm,s_red,s_copy,e_red,smatre,smatim,smat, &
               up,um,eigpre,eigpim,eigmre,eigmim,bpre,bpim,bmre,bmim)
!
    1050 format(t5,'----------------------------------------',/,&
        t7,'# target vectors:    ',i4,/,&
        t7,'# new vectors added: ',i4,/,&
        t7,'# converged vectors: ',i4,/,&
        t5,'----------------------------------------')
    return
end subroutine caslr_newcomplex_eff_driver
!
subroutine caslr_complex_eff_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                      apbmul_complex,ambmul_complex,spdmul_complex,smdmul_complex,lrprec_complex,eig,ceig,evec,ok)
!
!   main driver for the efficient solution to the following generalized 
!   eigenvalue problem:
!
!   / A   B  \ / Y  \     /  S   D  \ / Y  \ 
!   |        | |    | = w |         | |    |
!   \ B*  A* / \ Z* /     \ -D* -S* / \ Z* /
!
!   where A, B, S are symmetric matrices and D is antisymmetric.
!
!   If (w, Y, Z*) are a solution, then (-w, Z, Y*) is also a solution. 
!   Following J. Chem. Phys., 118, 522 (2003), we enforce this property in 
!   the iterative procedure by expanding the eigenvector as
!
!   (Y, Z*) = (b+, b+*) + (b-, -b-*)
!
!   This routine solves the associate problem 
!
!   /  S   D  \ / Y  \   1 / A   B  \ / Y  \ 
!   |         | |    | = - |        | |    |
!   \ -D* -S* / \ Z* /   w \ B*  A* / \ Z* /
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
!   TODO:  Here I need more subroutine to be implemented for complex lr 
!   Ab+, Bb+*, Ab-, -Bb-*, sigmab+, deltab+*, sigmab-, -deltab-*.
!   TODO: add the new subroutines to the call of caslr_complex ecc...
!
!   apbmul_complex:   external subroutine to compute (Ab+ + Bb+*) 
!
!   ambmul_complex:   external subroutine to compute (Ab- - Bb-*) 
!
!   spdmul_complex:   external subroutine to compute (Sb+ + Db+*) 
!
!   smdmul_complex:   external subroutine to compute (Sb- - Db-*) 
!
!   lrprec_complex:   external subroutine that applies a preconditioner.
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
    logical,                          intent(in)    :: verbose
    integer,                          intent(in)    :: n, n2, n_targ, n_max
    integer,                          intent(in)    :: max_iter, max_dav
    real(dp),                         intent(in)    :: tol
    real(dp),    dimension(n_max),    intent(inout) :: eig
    complex(dp), dimension(n_max),    intent(inout) :: ceig
    complex(dp), dimension(n2,n_max), intent(inout) :: evec
    logical,                          intent(inout) :: ok
    external                                        :: apbmul_complex, ambmul_complex, spdmul_complex, smdmul_complex, &
                                                       lrprec_complex
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
    complex(dp), allocatable :: vp(:,:), vm(:,:), lvp(:,:), lvm(:,:), bvp(:,:), bvm(:,:)
    complex(dp), allocatable :: rp(:,:), rm(:,:)
    real(dp),    allocatable :: rr(:,:), r_norm(:,:)
!
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
    real(dp),    allocatable :: up(:,:), um(:,:)
    complex(dp), allocatable :: cup(:,:), cum(:,:)
    complex(dp), allocatable :: eigp(:,:), eigm(:,:), bp(:,:), bm(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp),    allocatable :: s_red(:,:), s_copy(:,:), e_red(:)
    real(dp),    allocatable :: smat(:,:)
    complex(dp), allocatable :: csmat(:,:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2, dznrm2
    external              :: dcopy, dnrm2, dznrm2, dgemm, zgemm, dsyev, zheevd, zaxpy
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than max_dav = 10 is deemed acceptable.
!
    write(6,*) 'SMO-GD COMPLEX DIAGLIB '
!
    dim_dav = max(min_dav,max_dav)
    lda     = dim_dav*n_max
    lda2    = 2 * lda
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    allocate (work(lwork), cwork(lwork), tau(n_max), ctau(n_max), stat=istat)
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
    allocate (s_red(lda,lda), s_copy(lda,lda), e_red(lda2), smat(lda,lda), csmat(lda,lda), stat=istat)
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), cup(lda,n_max), cum(lda,n_max), eigp(n,n_max), eigm(n,n_max), &
              bp(n,n_max), bm(n,n_max), stat = istat)
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
    t_diag  = 0.0_dp
    t_ortho = 0.0_dp
    t_mv    = 0.0_dp
    t_tot   = 0.0_dp
    vp      = (0.0_dp,0.0_dp)
    vm      = (0.0_dp,0.0_dp)
!    up      = 0.0_dp
!    um      = 0.0_dp
!    cup     = (0.0_dp,0.0_dp)
!    cum     = (0.0_dp,0.0_dp)
    bvp     = (0.0_dp,0.0_dp)
    bvm     = (0.0_dp,0.0_dp)
    lvp     = (0.0_dp,0.0_dp)
    lvm     = (0.0_dp,0.0_dp)
!    smat    = 0.0_dp
!    csmat   = (0.0_dp,0.0_dp)
    ok      = .false.
    done    = .false.
!
    call get_time(t_tot)
!
!   move the guess into the expansion space.
!
    do i_eig = 1, n_max
      vp(:,i_eig) = evec(1:n,i_eig) + (evec(n+1:n2,i_eig))
      vm(:,i_eig) = evec(1:n,i_eig) - (evec(n+1:n2,i_eig))
    end do
!
!   orthogonalize the expansion space to the metric.
!
    call apbmul_complex(n,n_max,vp,lvp)
    call b_ortho_complex(n,n_max,vp,lvp)
    call ambmul_complex(n,n_max,vm,lvm)
!    write(6,*) 'vm_amb', vm
    call b_ortho_complex(n,n_max,vm,lvm)
!    write(6,*) 'Vp'
!    write(6,*) vp(:,2)
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
    1030 format(t5,'Complex SMOGD driver iterations (tol=',d10.2,'):',/, &
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
      call spdmul_complex(n,n_act,vp(1,i_beg),bvm(1,i_beg))
!      write(6,*) 'vm_prima', vm
      call smdmul_complex(n,n_act,vm(1,i_beg),bvp(1,i_beg))
!      write(6,*) 'Bvp', bvp
!      write(6,*) 'vm_dopo', vm
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call zgemm('c','n',ldu,ldu,n,cone,vm,n,bvm,n,czero,csmat,lda)
      smat = csmat + conjg(csmat)
!      smat = real(csmat) 
!      write(6,*) 'csmat', csmat + conjg(csmat)
!      write(6,*) 'smat', smat
      !write(6,*) 'smat'
      !call prtmat_complex(4,4,smat,3)!    3 is the free-format to correctly print a complex matrix
      !if (it .eq. 3) stop
      !write(6,*) 'S mat'
!      write(6,'(4f20.10)') transpose(smat(1:4,1:4))
      !write(6,*) transpose(smat(1:4,1:4))
!
!     save s, and assemble s^t s:
!
      s_red  = smat
      s_copy = 0.0_dp 
      call dgemm('t','n',ldu,ldu,ldu,one,s_red,lda,s_red,lda,zero,s_copy,lda)
!
!     diagonalize s^t s
!
      call get_time(t1)
      call dsyev('v','u',ldu,s_copy,lda,e_red,work,lwork,info)
      !print *, 'info:', info
      !print *, 'e_red:', e_red(1:lda)
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
      !write(6,*) 'eig', eig
!
!     assemble the symmetric and antisymmetric combinations (Y+Z) and (Y-Z)
!
      cup = up
      cum = um
      !write(6,*) 'cup', cup
      !write(6,*) 'up', up
      !write(6,*) 'cum', cum
      !write(6,*) 'um', um
      call zgemm('n','n',n,n_max,ldu,cone,vp,n,cup,lda,czero,eigp,n)
      call zgemm('n','n',n,n_max,ldu,cone,vm,n,cum,lda,czero,eigm,n)
!
!     assemble the current approximation to the eigenvectors
!
      do i_eig = 1, n_max
        evec(1:n,i_eig)    = eigp(:,i_eig) + eigm(:,i_eig)
        evec(n+1:n2,i_eig) = eigp(:,i_eig) - eigm(:,i_eig)
        evec(n+1:n2,i_eig) = conjg(evec(n+1:n2,i_eig))
      end do
!
!     compute the residuals, and their rms and sup norms:
!
      call zgemm('n','n',n,n_max,ldu,cone,bvp,n,cum,lda,czero,rp,n)
      call zgemm('n','n',n,n_max,ldu,cone,bvm,n,cup,lda,czero,rm,n)
      call zgemm('n','n',n,n_max,ldu,cone,lvp,n,cup,lda,czero,bp,n)
      call zgemm('n','n',n,n_max,ldu,cone,lvm,n,cum,lda,czero,bm,n)
!      
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        ceig(i_eig) = eig(i_eig)
        !call zaxpy(n,-dcmplx(eig(i_eig)),bp(:,i_eig),1,rp(:,i_eig),1)
        !call zaxpy(n,-dcmplx(eig(i_eig)),bm(:,i_eig),1,rm(:,i_eig),1)
        call zaxpy(n,-ceig(i_eig),bp(:,i_eig),1,rp(:,i_eig),1)
        call zaxpy(n,-ceig(i_eig),bm(:,i_eig),1,rm(:,i_eig),1)
        r_norm(1,i_eig) = (dznrm2(n,rp(:,i_eig),1) + dznrm2(n,rm(:,i_eig),1))/(ceig(i_eig)*sqrt(ctwo)*sqrtn)
        r_norm(2,i_eig) = (maxval(abs(rp(:,i_eig))) + maxval(abs(rm(:,i_eig))))/(sqrt(ctwo)*ceig(i_eig))
      end do
!      write(6,*) 'rp', rp
!      write(6,*) 'rm', rm
!      eig = real(ceig)
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
          ceig(i_eig) = cone/ceig(i_eig)
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
        call lrprec_complex(n,n_act,ceig(ind),rp(1,ind),rm(1,ind),vp(1,i_beg),vm(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        write(6,*) 'HERE pre b_ortho_vs_x 1'
        call b_ortho_vs_x_complex(n,ldu,n_act,vp,lvp,vp(1,i_beg))
        write(6,*) 'HERE post b_ortho_vs_x 1'
        call apbmul_complex(n,n_act,vp(1,i_beg),lvp(1,i_beg))
        call b_ortho_complex(n,n_act,vp(1,i_beg),lvp(1,i_beg))
        write(6,*) 'HERE pre b_ortho_vs_x 2'
        call b_ortho_vs_x_complex(n,ldu,n_act,vm,lvm,vm(1,i_beg))
        write(6,*) 'HERE post b_ortho_vs_x 2'
        call ambmul_complex(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call b_ortho_complex(n,n_act,vm(1,i_beg),lvm(1,i_beg))
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        write(6,*) 'HERE RESTART'
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
        vp = (0.0_dp,0.0_dp)
        vm = (0.0_dp,0.0_dp)
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max
          vp(:,i_eig) = evec(1:n,i_eig) + (evec(n+1:n2,i_eig))
          vm(:,i_eig) = evec(1:n,i_eig) - (evec(n+1:n2,i_eig))
        end do
!
        lvp = (0.0_dp,0.0_dp)
        lvm = (0.0_dp,0.0_dp)
!
        call apbmul_complex(n,n_max,vp,lvp)
        call b_ortho_complex(n,n_max,vp,lvp)
        call ambmul_complex(n,n_max,vm,lvm)
        call b_ortho_complex(n,n_max,vm,lvm)
!
        bvp   = (0.0_dp,0.0_dp)
        bvm   = (0.0_dp,0.0_dp)
        s_red = zero
        smat  = zero
        csmat = (0.0_dp,0.0_dp) 
!        
      end if
!      
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!      
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
    1000 format(t3,'timings for caslr_complex_eff (cpu/wall):   ',/, &
                t3,'  matrix-vector multiplications: ',2f12.4,/, &
                t3,'  diagonalization:               ',2f12.4,/, &
                t3,'  orthogonalization:             ',2f12.4,/, &
                t3,'                                 ',24('='),/,  &
                t3,'  total:                         ',2f12.4)
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
    deallocate(work,cwork,tau,ctau,vp,vm,lvp,lvm,bvp,bvm,rp,rm,rr,done,r_norm,s_red,s_copy,e_red,smat,up,um, & 
               cup,cum,eigp,eigm,bp,bm,csmat)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine caslr_complex_eff_driver
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
  subroutine davidson_newcomplex_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav, &
                                        shift,matvec,precnd,eig,evec,ok)
!
    logical,                      intent(in)    :: verbose
    integer,                      intent(in)    :: n, n_targ, n_max
    integer,                      intent(in)    :: max_iter, max_dav
    real(dp),                     intent(in)    :: tol, shift
    real(dp), dimension(n_max),   intent(inout) :: eig
    real(dp), dimension(n,n_max), intent(inout) :: evec !evecre, evecim
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
!    real(dp), allocatable :: spaceim(:,:), aspaceim(:,:), rim(:,:)
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
    integer               :: n2
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
!    n2 = 2*n
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
    space    = zero
    aspace   = zero
    a_red    = zero
    ok       = .false.
    done     = .false.
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess(n,n_max,evec)
!    call check_guess(n,n_max,evecim)
!
!   move the guess into the expansion space.
!
    call dcopy(n*n_max,evec,1,space,1)
!    call dcopy(n*n_max,evecim,1,spaceim,1)
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
!      call matvec(1,n,n_act,spaceim(1,i_beg+n_rst),aspaceim(1,i_beg+n_rst))
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
  end subroutine davidson_newcomplex_driver
!

  subroutine davidson_complex_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav, &
                             shift,matvec_complex,precnd_complex,eig,ceig,evec,ok)
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
    logical,                         intent(in)    :: verbose
    integer,                         intent(in)    :: n, n_targ, n_max
    integer,                         intent(in)    :: max_iter, max_dav
    real(dp),                        intent(in)    :: tol, shift
    real(dp),    dimension(n_max),   intent(inout) :: eig
    complex(dp), dimension(n_max),   intent(inout) :: ceig
    complex(dp), dimension(n,n_max), intent(inout) :: evec
    logical,                         intent(inout) :: ok
    external                                       :: matvec_complex, precnd_complex
!
!   lapack: zheevd 
!
    integer                  :: lrwork, liwork
    integer,     allocatable :: iwork(:) 
    real(dp),    allocatable :: w(:), rwork(:)
!    complex(dp), allocatable :: work(:)
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
    complex(dp)           :: xx(1)
!
!   arrays to control convergence and orthogonalization
!
    logical,  allocatable :: done(:)
!
!   expansion spaces, residuals and their norms.
!
    complex(dp), allocatable :: space(:,:), aspace(:,:), r(:,:)
    real(dp),    allocatable :: r_norm(:,:)
!
!   subspace matrix and eigenvalues.
!
    complex(dp), allocatable :: a_red(:,:), a_copy(:,:)
    real(dp),    allocatable :: e_red(:)
!
!   restarting variables
!
    integer               :: n_rst
    logical               :: restart
!
!   external functions:
!   ===================
!
    real(dp)              :: dnrm2, dznrm2
    external              :: dcopy, dznrm2, dnrm2, dgemm, dsyev, zgemm, zheevd, zcopy 
!
    write(6,*) 'DAVIDSON_COMPLEX_DRIVER RE DEI POPOLI OPPRESSI'
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
    !write(6,*) 'lwork at the beginning:', lwork
    allocate (cwork(lwork), rwork(1), iwork(1), tau(n_max), ctau(n_max), stat=istat)
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
    space    = (0.0_dp,0.0_dp)
    aspace   = (0.0_dp,0.0_dp)
    a_red    = (0.0_dp,0.0_dp)
    ok       = .false.
    done     = .false.
!
    call get_time(t_tot)
!
!   check whether we have a guess for the eigenvectors in evec, and
!   whether it is orthonormal.
!   if evec is zero, create a random guess.
!
    call check_guess_complex(n,n_max,evec)
!
!   move the guess into the expansion space.
!
    call zcopy(n*n_max,evec,1,space,1)
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
      call matvec_complex(n,n_act,space(1,i_beg+n_rst),aspace(1,i_beg+n_rst))
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call zgemm('c','n',ldu,n_act,n,cone,space,n,aspace(1,i_beg+n_rst),n,czero,a_red(1,i_beg+n_rst),lda)
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
      call zheevd('v','u',ldu,a_copy,lda,e_red,cwork,-1,rwork,-1,iwork,-1,info)
      write(6,*) 'info1', info
      lwork = int(cwork(1))
      lrwork = int(rwork(1))
      liwork = int(iwork(1))
      deallocate(cwork, rwork, iwork)
      allocate (cwork(lwork), rwork(lrwork), iwork(liwork))
      call get_time(t1)
      !write(6,*) 'a red'
      !write(6,*) a_copy(1:10,1:10)
      call zheevd('v','u',ldu,a_copy,lda,e_red,cwork,lwork,rwork,lrwork,iwork,liwork,info)
      write(6,*) 'info reduced matrix', info
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      eig = e_red(1:n_max)
      !write(6,*) 'e_red', e_red
!
      call zgemm('n','n',n,n_max,ldu,cone,space,n,a_copy,lda,czero,evec,n)
      !write(6,*) 'evec', imag(evec) 
!
!     compute the residuals, and their rms and sup norms:
!
      call zgemm('n','n',n,n_max,ldu,cone,aspace,n,a_copy,lda,czero,r,n)
!
      do i_eig = 1, n_targ
!
!       if the eigenvalue is already converged, skip it.
!
        if (done(i_eig)) cycle
!
        ceig(i_eig) = eig(i_eig)
        call zaxpy(n,-ceig(i_eig),evec(:,i_eig),1,r(:,i_eig),1)
        r_norm(1,i_eig) = dznrm2(n,r(:,i_eig),1)/sqrtn
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
        call precnd_complex(n,n_act,-ceig(ind),r(1,ind),space(1,i_beg))
!
!       orthogonalize the new vectors to the existing ones and then
!       orthonormalize them.
!
        call get_time(t1)
        call ortho_vs_x_complex(.false.,n,ldu,n_act,space,space(1,i_beg),xx,xx,0)
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
      else
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
        n_act = n_max
        space = (0.0_dp,0.0_dp)
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        call zcopy(n_max*n,evec,1,space,1)
        aspace = (0.0_dp,0.0_dp)
        a_red  = (0.0_dp,0.0_dp)
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
    deallocate (cwork, rwork, iwork, ctau, tau, space, aspace, r, done, r_norm, a_red, a_copy, e_red)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
    return
  end subroutine davidson_complex_driver
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
  subroutine ortho_complex(do_other,n,m,u,w)
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
    logical,                      intent(in)    :: do_other
    integer,                      intent(in)    :: n, m
    complex(dp),  dimension(n,m), intent(inout) :: u, w
!
!   local scratch
!   =============
!
    complex(dp), allocatable :: v(:,:)
!
!   external functions:
!   ===================
!
    external zgeqrf, ztrsm
!
    allocate (v(n,m))
    v = u
!    write(6,*) 'lwork', lwork
    call zgeqrf(n,m,u,n,ctau,cwork,lwork,info)
    write(6,*) 'info zgeqrf', info
    !here
!
    call ztrsm('r','u','n','n',n,m,cone,u,n,v,n)
!
    if (do_other) call ztrsm('r','u','n','n',n,m,cone,u,n,w,n)
!
    u = v
!
    deallocate (v)
    return
  end subroutine ortho_complex
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
  subroutine b_ortho_complex(n,m,u,bu)
    implicit none
!
!   b-orthogonalize m complex vectors of lenght n using the cholesky decomposition
!   of the overlap matrix.
!   this is in principle not a good idea, as the u'bu matrix can be very
!   ill-conditioned, independent of how bas is b, and only works if x is
!   already orthonormal. 
!
!   arguments:
!   ==========
!
    integer,                        intent(in)    :: n, m
    complex(dp),  dimension(n,m),   intent(inout) :: u, bu
!
!   local variables
!   ===============
!
    integer                  :: info, i, j
    complex(dp), allocatable :: metric(:,:), csigma(:), u_svd(:,:), vt_svd(:,:), &
                                temp(:,:)
    real(dp),    parameter   :: tol_svd = 1.0e-5_dp
    real(dp),    allocatable :: rwork(:), sigma(:), metric_re(:,:)
    logical,     parameter   :: use_svd = .false.
!    intrinsic                :: zsqrt
!
!   external functions:
!   ===================
!
    external dpotrf, zpotrf, ztrsm, zgemm, zgesvd
!
    allocate (metric(m,m), metric_re(m,m))
!
    print *, 'metric pre'
    call prtmat_complex(4,n,metric,3)!    3 is the free-format to correctly print a complex matrix
    call zgemm('c','n',m,m,n,cone,u,n,bu,n,czero,metric,m)
    !metric = metric  + conjg(metric
    metric = 2.0_dp*real(metric)
    write(6,*) 'metric post', metric
    !print *, 'metric post'
    !call prtmat(4,n,metric_re,2)
!
    if (use_svd) then
!
!     debug option: use svd to b-orthonormalize, by computing
!     b**(-1/2)
!
      allocate (sigma(m), csigma(m), u_svd(m,m), vt_svd(m,m), temp(n,m), rwork(5*m))
      call zgesvd('a','a',m,m,metric,m,sigma,u_svd,m,vt_svd,m,cwork,lwork,rwork,info)
      write(6,*) 'b_ortho_complex HEREEEEE!!!!!!', info
      if (info .ne. 0) write(6,*) 'INFO ZGESVD IS NOT ZERO, b_ortho_complex', info
!      sigma = real(csigma)
!
!     compute sigma**(-1/2)
!
      do i = 1, m
        if (sigma(i) .gt. tol_svd) then
          !sigma(i) = dcmplx(1.0_dp,1.0_dp)/sqrt(sigma(i))
          sigma(i) = 1.0_dp/sqrt(sigma(i))
        else
          sigma(i) = 0.0_dp
        end if
      end do    
!
!     compute metric ** (-1/2). first, compute sigma ** (-1/2) vt
!
      metric = (0.0_dp,0.0_dp)
      do i = 1, m
        do j = 1, m
          metric(j,i) = metric(j,i) + sigma(j)*vt_svd(j,i)
        end do
      end do
!
!     now, multiply for u:
!
      vt_svd = metric
      call zgemm('n','n',m,m,m,cone,u_svd,m,vt_svd,m,czero,metric,m)
!
!     metric contains s ** (-1/2), and projects out directions corresponding
!     to pathological singular values. 
!     orthogonalize u and bu:
!
      call zgemm('n','n',n,m,m,cone,u,n,metric,m,czero,temp,n)
      u = temp
      call zgemm('n','n',n,m,m,cone,bu,n,metric,m,czero,temp,n)
      bu = temp
!
      deallocate (sigma, csigma, u_svd, vt_svd, temp)
    else
!
!     compute the cholesky factorization of the metric.
!
!      call dpotrf('l',m,metric_re,m,info)
!      metric = (0.0_dp,0.0_dp)
!      metric = metric_re
      call zpotrf('l',m,metric,m,info)
!
!     get u * l^-T and bu * l^-T
!
      call ztrsm('r','l','c','n',n,m,cone,metric,m,u,n)
      call ztrsm('r','l','c','n',n,m,cone,metric,m,bu,n)
!
    end if
!
    deallocate (metric, metric_re)
!
    return
  end subroutine b_ortho_complex
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
  subroutine ortho_cd_newcomplex(do_other,n,m,u,w,ok)
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
      call dgemm('t','n',m,m,n,two,u,n,u,n,zero,metric,m)
      !call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
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
      !call dtrsm('r','l','t','n',n,m,two,metric,m,u,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,u,n)
      if (do_other) call dtrsm('r','l','t','n',n,m,one,metric,m,w,n)
!
!     check that the vectors in v are really orthogonal
!
      call dgemm('t','n',m,m,n,two,u,n,u,n,zero,metric,m)
      !call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
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
  end subroutine ortho_cd_newcomplex
!
  subroutine ortho_cd_complex(do_other,n,m,u,w,ok)
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
    logical,                      intent(in)    :: do_other
    integer,                      intent(in)    :: n, m
    complex(dp),  dimension(n,m), intent(inout) :: u, w
    logical,                      intent(inout) :: ok
!
!   local variables
!   ===============
!
    integer               :: it, it_micro, i
    real(dp)              :: metric_norm, dnrm2, dznrm2, alpha, unorm 
    real(dp)              :: shift !dc
    complex(dp)           :: calpha!, shift
    logical               :: macro_done, micro_done
    integer, parameter    :: maxit = 10

!   local scratch
!   =============
!
    complex(dp), allocatable :: metric(:,:), msave(:,:)
    real(dp), allocatable :: ureal(:,:)
!
!   external functions:
!   ===================
!
    external              :: zgemm, zgeqrf, ztrsm, dnrm2, dznrm2, zpotrf
!
!   get memory for the metric.
!
    allocate (metric(m,m), msave(m,m), ureal(m,m))
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
      call zgemm('c','n',m,m,n,cone,u,n,u,n,czero,metric,m)
      !ceci
      !metric = metric + conjg(metric)
      !ureal = u
      !print *, 'v+ * v+:'
      !print '(10f15.8)', transpose(ureal)
!      write(6,*) 'metric'
!      write(6,*) metric
      msave = metric
!
!   compute the cholesky factorization of the metric.
!
      call zpotrf('l',m,metric,m,info)
!      if (info .ne. 0) then
!        write(6,*) 'INFO ZPOTRF IS NOT ZERO, ortho_cd_complex', info
!      end if
!
!     if zpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!
      if (info.ne.0) then
!
        calpha     = (100.0_dp,0.0_dp)
        alpha      = 100.0_dp
        unorm      = dznrm2(n*m,u,1)
        it_micro   = 0
        micro_done = .false.
!
!       add larger and larger shifts to the diagonal until zpotrf manages to factorize it.
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
          shift = max(abs(epsilon(one)*alpha*unorm),tol_ortho)
          metric = msave
          call diag_shift_complex(m,shift,metric)
          call zpotrf('l',m,metric,m,info)
          !calpha = calpha * dcmplx(10.0_dp,10.0_dp)
          calpha = calpha * (10.0_dp,0.0_dp)
          alpha  = alpha * 10.0_dp
          micro_done = info.eq.0
        end do
!
      end if
!
!     orthogonalize u by computing a solution to u(ortho) l^t = u
!     if required, apply the same transformation to w.
!    
      call ztrsm('r','l','c','n',n,m,cone,metric,m,u,n)
      if (do_other) call ztrsm('r','l','c','n',n,m,cone,metric,m,w,n)
!
!     check that the vectors in v are really orthogonal
!
      call zgemm('c','n',m,m,n,cone,u,n,u,n,czero,metric,m)
      !ceci
      !metric = metric + conjg(metric)
!      print *, metric
      metric_norm = abs(dznrm2(m*m,metric,1) - sqrt(dble(m)))
      macro_done  = metric_norm .lt. tol_ortho
!      write(6,*) 'metric_norm - tol_ortho', metric_norm, tol_ortho
    end do
!
    100 format(t3,'ortho_cd_complex failed with the following error:',a)
!
    ok = .true.
    !do i = 1, m
    !  u(:,i) = u(:,i) / dznrm2(n,u(:,i),1)
    !end do
!
    deallocate (metric)
    return
  end subroutine ortho_cd_complex
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
    if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
    ! I put .and. instead of .or.
    if (.not. ok .and. useqr) call ortho(do_other,n,k,u,au)
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
      ! I put .and. instead of .or.
      if (.not. ok .and. useqr) call ortho(do_other,n,k,u,au)
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
  subroutine ortho_vs_x_newcomplex(do_other,n,m,k,x,u,ax,au)
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
    if (.not. useqr) call ortho_cd_newcomplex(do_other,n,k,u,au,ok)
    ! I put .and. instead of .or.
    if (.not. ok .and. useqr) call ortho(do_other,n,k,u,au)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^t u)
!
      call dgemm('t','n',m,k,n,two,x,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
      if (do_other) call dgemm('n','n',n,k,m,-one,ax,n,xu,m,one,au,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd_newcomplex(do_other,n,k,u,au,ok)
      ! I put .and. instead of .or.
      if (.not. ok .and. useqr) call ortho(do_other,n,k,u,au)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      call dgemm('t','n',m,k,n,two,x,n,u,n,zero,xu,m)
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
  end subroutine ortho_vs_x_newcomplex
!
  subroutine ortho_vs_x_4n(do_other,n,m,k,x,u,ax,au,iwhat)
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
    integer,                   intent(in)    :: n, m, k, iwhat
    real(dp),  dimension(n,m), intent(inout) :: x, ax
    real(dp),  dimension(n,k), intent(inout) :: u, au
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it, i_eig, i
    real(dp)               :: xu_norm
    real(dp),  allocatable :: xu(:,:), up(:), down(:)
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
    allocate (xu(m,k), up(n/4), down(n/4))
    done = .false.
    it   = 0
!check    if (iwhat .eq. 1) then
!check      write(6,*) 'before ortho'
!check      call check_sym_4n(1,n,m,x)
!check    else if (iwhat .eq. 2) then
!check      write(6,*) 'before ortho'
!check      call check_sym_4n(2,n,m,x)
!check    end if
    !do i_eig = 1, m
    !  print *, 'x', i_eig
    !  do i = 1, n/4
    !    print '(4f15.8)', x(i,i_eig), x(i+n/4,i_eig), x(n/2+i,i_eig), x(3*n/4+i,i_eig)
    !  end do
    !end do
    !do i_eig = 1, k
    !  print *, 'u', i_eig
    !  do i = 1, n/4
    !    print '(4f15.8)', u(i,i_eig), u(i+n/4,i_eig), u(n/2+i,i_eig), u(3*n/4+i,i_eig)
    !  end do
    !end do
!
!   start with an initial orthogonalization to improve conditioning.
!
    if (iwhat .eq. 1) then
!check      write(6,*) 'U before ortho'
!check      call check_sym_4n(1,n,k,u)
      do i = 1, k
        up   = (u(1:n/4,i) + u(n/4+1:n/2,i))/2.0d0
        down = (u(n/2+1:3*n/4,i) - u(3*n/4+1:n,i))/2.0d0
        u(1:n/4,i)     = up
        u(n/4+1:n/2,i) = up
        u(n/2+1:3*n/4,i) = down
        u(3*n/4+1:n,i)   = -down
      end do
!check      write(6,*) 'U before ortho cleaned'
!check      call check_sym_4n(1,n,k,u)
    else if (iwhat .eq. 2) then
!check      write(6,*) 'U before ortho'
!check      call check_sym_4n(2,n,k,u)
      do i = 1, k
        up   = (u(1:n/4,i) - u(n/4+1:n/2,i))/2.0d0
        down = (u(n/2+1:3*n/4,i) + u(3*n/4+1:n,i))/2.0d0
        u(1:n/4,i)     = up
        u(n/4+1:n/2,i) = -up
        u(n/2+1:3*n/4,i) = down
        u(3*n/4+1:n,i)   = down
      end do
!check      write(6,*) 'U before ortho cleaned'
!check      call check_sym_4n(2,n,k,u)
    end if
    if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
    ! I put .and. instead of .or.
    if (.not. ok .and. useqr) call ortho(do_other,n,k,u,au)
!check    if (iwhat .eq. 1) then
!check      write(6,*) 'U after ortho'
!check      call check_sym_4n(1,n,k,u)
!check    else if (iwhat .eq. 2) then
!check      write(6,*) 'U after ortho'
!check      call check_sym_4n(2,n,k,u)
!check    end if
    !if (.not. ok .or. useqr) call ortho(do_other,n,k,u,au)
    !do i_eig = 1, k
    !  print *, 'u post ortho', i_eig
    !  do i = 1, n/4
    !    print '(4f15.8)', u(i,i_eig), u(i+n/4,i_eig), u(n/2+i,i_eig), u(3*n/4+i,i_eig)
    !  end do
    !end do
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
      !do i_eig = 1, k
      !  print *, 'u post GS', i_eig
      !  do i = 1, n/4
      !    print '(4f15.8)', u(i,i_eig), u(i+n/4,i_eig), u(n/2+i,i_eig), u(3*n/4+i,i_eig)
      !  end do
      !end do
!
      if (do_other) call dgemm('n','n',n,k,m,-one,ax,n,xu,m,one,au,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd(do_other,n,k,u,au,ok)
      ! I put .and. instead of .or.
      if (.not. ok .and. useqr) call ortho(do_other,n,k,u,au)
      !if (.not. ok .or. useqr) call ortho(do_other,n,k,u,au)
      !do i_eig = 1, k
      !  print *, 'u post ortho GS', i_eig
      !  do i = 1, n/4
      !    print '(4f15.8)', u(i,i_eig), u(i+n/4,i_eig), u(n/2+i,i_eig), u(3*n/4+i,i_eig)
      !  end do
      !end do
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
    deallocate(xu, up, down)
!
    return
  end subroutine ortho_vs_x_4n
!
  subroutine ortho_vs_x_complex(do_other,n,m,k,x,u,ax,au,bsign)
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
    logical,                      intent(in)    :: do_other
    integer,                      intent(in)    :: n, m, k, bsign
    complex(dp),  dimension(n,m), intent(in)    :: x, ax
    complex(dp),  dimension(n,k), intent(inout) :: u, au
!
!   local variables:
!   ================
!
    logical                   :: done, ok
    integer                   :: it, i, j
    real(dp)                  :: xu_norm
    complex(dp), allocatable  :: xu(:,:)
!
!   external functions:
!   ===================
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2, dznrm2
    external               :: dnrm2, dgemm, dznrm2, zgemm
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
    if (.not. useqr) call ortho_cd_complex(do_other,n,k,u,au,ok)
    if (.not. ok .and. useqr) call ortho_complex(do_other,n,k,u,au)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^H u)
!
      call zgemm('c','n',m,k,n,cone,x,n,u,n,czero,xu,m)
      !ceci
!      write(6,*) 'bsign', bsign
      if (bsign == 0) then
        xu = xu 
      else if (bsign == 1) then
        xu = xu - conjg(xu)
      else if (bsign == 2) then
        do i = 1, m/2
          xu(i,:) = xu(i,:) + conjg(xu(i,:))
          xu(i+m/2,:) = xu(i+m/2,:) - conjg(xu(i+m/2,:))
        end do
      else if (bsign == 3) then 
        do i = 1, (m-k)/2
          xu(i,:) = xu(i,:) - conjg(xu(i,:))
          xu(i+(m-k)/2,:) = xu(i+(m-k)/2,:) + conjg(xu(i+(m-k)/2,:))
        end do
        do i = m-k+1, m
          xu(i,:) = xu(i,:) - conjg(xu(i,:))
        end do
      else if (bsign == 4) then
        xu = xu + conjg(xu) 
      else
        write(6,*) 'invalid bsign entry, aborting'
      end if
!      if (bsign) xu = xu + conjg(xu)
!      if (.not. bsign) xu = xu - conjg(xu)
      call zgemm('n','n',n,k,m,-cone,x,n,xu,m,cone,u,n)
!
      if (do_other) call zgemm('n','n',n,k,m,-cone,ax,n,xu,m,cone,au,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd_complex(do_other,n,k,u,au,ok)
      if (.not. ok .and. useqr) call ortho_complex(do_other,n,k,u,au)
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      call zgemm('c','n',m,k,n,cone,x,n,u,n,czero,xu,m)
      !ceci
      if (bsign == 0) then
        xu = xu  
      else if (bsign == 1) then
        xu = xu - conjg(xu)   
      else if (bsign == 2) then
        do i = 1, m/2
          xu(i,:) = xu(i,:) + conjg(xu(i,:))
          xu(i+m/2,:) = xu(i+m/2,:) - conjg(xu(i+m/2,:))
        end do
      else if (bsign == 3) then 
        do i = 1, (m-k)/2
          xu(i,:) = xu(i,:) - conjg(xu(i,:))
          xu(i+(m-k)/2,:) = xu(i+(m-k)/2,:) + conjg(xu(i+(m-k)/2,:))
        end do
        do i = m-k+1, m
          xu(i,:) = xu(i,:) - conjg(xu(i,:))
        end do
      else if (bsign == 4) then
        xu = xu + conjg(xu) 
      else 
        write(6,*) 'invalid bsign entry, aborting'
      end if
      !tnwrite(6,*) "Xu values:"
      !tndo i = 1, k
      !tn  do j = 1, m
      !tn    write(6,*) xu(j,i)
      !tn  end do
      !tnend do
!
!      if (bsign) xu = xu + conjg(xu)
!      if (.not. bsign) xu = xu - conjg(xu)
      xu_norm = dznrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of ortho_vs_x_complex'
    end do
!
    deallocate(xu)
    !do i = 1, k
    !  u(:,i) = u(:,i) / dznrm2(n,u(:,i),1)
    !end do
!
    return
  end subroutine ortho_vs_x_complex
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
  subroutine b_ortho_vs_x_complex(n,m,k,x,bx,u)
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
    integer,                      intent(in)    :: n, m, k
    complex(dp),  dimension(n,m), intent(in)    :: x, bx
    complex(dp),  dimension(n,k), intent(inout) :: u
!
!   local variables:
!   ================
!
    logical                   :: done, ok
    integer                   :: it
    real(dp)                  :: xu_norm
    complex(dp)               :: xx(1)
    complex(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   ===================
!
    intrinsic              :: random_number
    real(dp)               :: dnrm2, dznrm2
    external               :: dnrm2, zgemm, dznrm2
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
    if (.not. useqr) call ortho_cd_complex(.false.,n,k,u,xx,ok)
    if (.not. ok .or. useqr) call ortho_complex(.false.,n,k,u,xx)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (bx^t u)
!
      call zgemm('c','n',m,k,n,cone,bx,n,u,n,czero,xu,m)
      call zgemm('n','n',n,k,m,-cone,x,n,xu,m,cone,u,n)
!      call zgemm('c','n',m,k,n,cone,u,n,bx,n,czero,xu,m)
!      call zgemm('n','n',n,k,m,-cone,x,n,xu,m,cone,u,n)
!
!     now, orthonormalize u.
!
      !write(6,*) 'HEREEEEEEEEE vs x1'
      if (.not. useqr) call ortho_cd_complex(.false.,n,k,u,xx,ok)
      !stop
      !write(6,*) 'HEREEEEEEEEE vs x2'
      if (.not. ok .or. useqr) call ortho_complex(.false.,n,k,u,xx)
      !write(6,*) 'HEREEEEEEEEE vs x3'
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      !write(6,*) 'HEREEEEEEEEE vs pre zgemm'
      call zgemm('c','n',m,k,n,cone,bx,n,u,n,czero,xu,m)
      !call zgemm('c','n',m,k,n,cone,u,n,bx,n,czero,xu,m)
      !write(6,*) 'HEREEEEEEEEE vs post zgemm'
      xu_norm = dznrm2(m*k,xu,1)
      write(6,*) 'xu - ortho - it', xu_norm, tol_ortho, it
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
     if (it.gt.maxit) stop ' catastrophic failure of b_ortho_vs_x_complex'
     ! if (it.gt.maxit) write(6,*) 'WARNING: b-orthogonalization converged in more than 10 iterations' 
    end do
!
    deallocate(xu)
!
    return
  end subroutine b_ortho_vs_x_complex
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
  subroutine diag_shift_complex(n,shift,a)
    implicit none
!  
!   add shift to the diagonal elements of the metric a
!  
    integer,                     intent(in)    :: n
    real(dp),                    intent(in)    :: shift
!    complex(dp),                 intent(in)    :: shift
    complex(dp), dimension(n,n), intent(inout) :: a
!  
    integer :: i
!  
    do i = 1, n
      a(i,i) = a(i,i) + shift
    end do
!  
    return
  end subroutine diag_shift_complex
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
  subroutine check_guess_complex(n,m,evec)
    implicit none
    integer,                     intent(in)    :: n, m
    complex(dp), dimension(n,m), intent(inout) :: evec
!
    integer                   :: i, j, istat
    real(dp)                  :: fac, diag_norm, out_norm
    complex(dp)               :: xx(1)
!
    complex(dp), allocatable  :: coverlap(:,:)
    real(dp),    allocatable  :: overlap(:,:)
    real(dp)                  :: dnrm2, dznrm2
    external                  :: dnrm2, dznrm2
!    
    integer,     dimension(4) :: vec
    real(dp),    dimension(4) :: vecr
!
!   check whether evec is zero.
!
    fac = dznrm2(n*m,evec,1)
    if (fac.eq.zero) then
!
!     no luck. make a random guess, then orthonormalize it.
!
      vecr = 0.0_dp
      call random_number(vecr)
      vec = floor(101*vecr,dp) 
      if (mod(vec(4),2) .ne. 1) vec(4)=vec(4)+1
      call zlarnv(4,vec,n*m,evec)  
      call ortho_complex(.false.,n,m,evec,xx)
    else
!
!     compute the overlap and check that the vectors are orthonormal.
!
      allocate (overlap(m,m), coverlap(m,m), stat=istat)
      call check_mem(istat)
      call zgemm('c','n',m,m,n,cone,evec,n,evec,n,czero,coverlap,m)
!      overlap = real(coverlap)
      diag_norm = zero
      out_norm  = zero
      do i = 1, m
!        diag_norm = diag_norm + dznrm2(1,coverlap(i,i),1)**2
        diag_norm = diag_norm + conjg(coverlap(i,i))*coverlap(i,i)
        do j = 1, i-1
          out_norm = out_norm + conjg(coverlap(i,j))*coverlap(i,j)
        end do
      end do
!
      diag_norm = diag_norm/real(m,dp)
!
      if (diag_norm .ne. one .or. out_norm.ne.zero) then
!
!       orthogonalize the guess:
!
        call ortho_complex(.false.,n,m,evec,xx)
      end if
!
      deallocate (overlap, coverlap, stat = istat)
      call check_mem(istat)
    end if
!
    return
  end subroutine check_guess_complex
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
  subroutine prtmat_complex(n,lda,mat,iform)
    integer,                         intent(in) :: n, lda, iform
    complex(dp), dimension(lda,lda), intent(in) :: mat
!
    integer :: i
!
    100 format(t3,20f12.6)
    200 format(t3,20d12.4)
!
    write(6,*) 'REAL PART: '
    do i = 1, n
      if (iform.eq.1) write(6,100) mat(i,1:n)
      if (iform.eq.2) write(6,200) mat(i,1:n)
      if (iform.eq.3) write(6,*) dble(mat(i,1:n))
    end do
    if (iform.eq.3) write(6,*) 
    write(6,*) 'IMAGINARY PART: '
    do i = 1, n
      if (iform.eq.3) write(6,*) dimag(mat(i,1:n))
    end do
    return
  end subroutine prtmat_complex
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
!
  subroutine check_sym_4n(iwhat,n,m,v)
    implicit none 
    integer,                  intent(in)  :: iwhat, n, m
    real(dp), dimension(n,m), intent(in)  :: v
    integer                               :: i
!
!   here we want to check the orthogonality of the expansion subspace vectors V(+) or V(-)
!
    if (iwhat .eq. 1) then
!
!     check the V(+)
!   
      write(6,*) 'CHECKING V(+)'
      do i = 1, m
        write(6,*)           'column number:      ', i 
        write(6,'(a,d20.6)') 'norm(br+ - br+):    ', norm2(v(1:n/4,i) - v(n/4+1:n/2,i))
        write(6,'(a,d20.6)') 'norm(bi+ + (-bi+)): ', norm2(v(n/2+1:3*n/4,i) + v(3*n/4+1:n,i))
      end do 
    else if (iwhat .eq. 2) then 
!
!     check the V(-)
!   
      write(6,*) 'CHECKING V(-)'
      do i = 1, m
        write(6,*)           'column number:      ', i 
        write(6,'(a,d20.6)') 'norm(br- + (-br-)): ', norm2(v(1:n/4,i) + v(n/4+1:n/2,i))
        write(6,'(a,d20.6)') 'norm(bi- - bi-):    ', norm2(v(n/2+1:3*n/4,i) - v(3*n/4+1:n,i))
      end do 
    else
      write(6,*) 'invalid entry for subroutine check_sym_4n'
    end if
  end subroutine check_sym_4n 
!
end module diaglib
