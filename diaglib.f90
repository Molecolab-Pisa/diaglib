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
! FOR EACH ORTHOGONALIZATION SUBROUTINE EXISTS A "COMPLEX" COUNTERART WHICH HAS 
! A DIFFERNT INTERFACE, SO THAT WE PASS BOTH THE REAL AND IMAGINARY PART OF 
! A CERTAIN OBJECT. HOW THE SUBROUTINE WORKS DOES NOT CHANGE. 
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
  real(dp),    parameter    :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10.0_dp
!
! convergence thresholds for orthogonalization
!
  real(dp), parameter       :: tol_ortho = 1.0e-12_dp
! 
! memory and info for lapack routines
!
  integer                   :: lwork, info
  real(dp),    allocatable  :: work(:), tau(:), rwork(:)
!
! timings:
!
  real(dp)                  :: t1(2), t2(2), t_diag(2), t_ortho(2), &
                               t_mv(2), t_tot(2)
!
! subroutines:
! ============
!
  public :: lobpcg_driver, davidson_driver, gen_david_driver, caslr_driver, caslr_eff_driver, caslr_complex_driver, & 
            caslr_complex_eff_driver, ortho, b_ortho, ortho_cd, ortho_vs_x, b_ortho_vs_x, prtmat, check_sym_4n, & 
            ortho_complex, ortho_cd_complex, b_ortho_vs_x_complex, b_ortho_complex
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
  subroutine caslr_complex_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                                  apbmul,ambmul,spdmul,smdmul,lrprec,eig,evecre,evecim,ok)
!
! The following driver is based on the traditional algorithm (Olsen's), 
! which solves the linear response equations defined as:
!
!   / A   B \ / Y \     /  S   D \ / Y \ 
!   |       | |   | = w |        | |   |
!   \ B*  A*/ \ Z*/     \ -D* -S*/ \ Z*/
!
! where A, B, S are hermitian matrices and D is antisymmetric. 
! Here, we want to solve the complex linear response equations
! using a real algebra formalism.
! The main idea is to apply a sym-antisym separation to each object
! (matrices and vectors) in order to remove the complex algebra, 
! and therefore obtaining a 4nx4n matrix-problem: 
! 
! /  Ar  Br -Ai -Bi \ /  Yr \     /  Sr  Dr -Si -Di \ /  Yr \
! |  Br  Ar  Bi  Ai | |  Zr | = w | -Dr -Sr -Di -Si | |  Zr |
! |  Ai  Bi  Ar  Br | |  Yi |     |  Si  Di  Sr  Dr | |  Yi |
! \ -Bi -Ai  Br  Ar / \ -Zi /     \  Di  Si -Dr -Sr / \ -Zi /
!
! As a consequence, the subspace vectors are so defined: 
! 
! b(+) \in R^4n = (b_r+  b_r+ b_i+ -b_i+),
! b(-) \in R^4n = (b_r- -b_r- b_i-  b_i-)
!
! This formalism ensures the orthogonality of the expansion vectors.
!
! Solving the augmented form of the LR-eqs leads to degenerate coupled-eigenvalues, 
! which belong to different eigenvector: 
!
! X1 = ( Yr Zr Yi -Zi), 
! X2 = (-Yi Zi Yr  Zr).
!
! As was already done in the "REAL" solver (caslr_driver) we only manage to deal with
! a reduced form of the expansion vectors; 
! hence, we use only b_r+/b_i+, b_r-/b_i- to build the expansion subspace.
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
    external                                     :: apbmul, ambmul, spdmul, smdmul, lrprec
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
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
!   real
    real(dp), allocatable :: up(:,:), um(:,:), eigpre(:,:), eigmre(:,:), bpre(:,:), bmre(:,:)
!   imaginary
    real(dp), allocatable :: eigpim(:,:), eigmim(:,:), bpim(:,:), bmim(:,:)
!
!   subspace matrix and eigenvalues.
!
    real(dp), allocatable :: a_red(:,:), a_copy(:,:), s_red(:,:), s_copy(:,:), e_red(:)
    real(dp), allocatable :: epmat(:,:),   emmat(:,:),   smat(:,:)
    real(dp), allocatable :: epmatre(:,:), emmatre(:,:), smatre(:,:)
    real(dp), allocatable :: epmatim(:,:), emmatim(:,:), smatim(:,:)
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
              emmatre(lda,lda), emmatim(lda,lda), smat(lda,lda), smatre(lda,lda), smatim(lda,lda), stat=istat) 
    call check_mem(istat)
!
!   allocate memory for the plus and minus eigenvector components:
!
    allocate (up(lda,n_max), um(lda,n_max), eigpre(n,n_max), eigpim(n,n_max), eigmre(n,n_max), eigmim(n,n_max), & 
              bpre(n,n_max), bpim(n,n_max), bmre(n,n_max), bmim(n,n_max), stat = istat)
    call check_mem(istat)
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
    a_red     = zero
    ok        = .false.
    done      = .false.
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
      vpre(:,i_eig+1) = (evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vmre(:,i_eig+1) = (evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!     imaginary 
      vpim(:,i_eig)   = ( evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig)   = ( evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vpim(:,i_eig+1) = (-evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig+1) = (-evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
    end do
!
    call ortho_cd_complex(.false.,n,n_max,vpre,vpim,xx,ok)
    call ortho_cd_complex(.false.,n,n_max,vmre,vmim,xx,ok)
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
      call dgemm('t','n',ldu,ldu,n,one,vpre,n,lvpre,n,zero,epmatre,lda)   !Ere+=(Vre+)^T*sigma_re+
      call dgemm('t','n',ldu,ldu,n,one,vpim,n,lvpim,n,zero,epmatim,lda)   !Eim+=(Vim+)^T*sigma_im+
      epmat = epmatre + epmatim
!      
      call dgemm('t','n',ldu,ldu,n,one,vmre,n,lvmre,n,zero,emmatre,lda)   !Ere-=(Vre-)^T*sigma_re-
      call dgemm('t','n',ldu,ldu,n,one,vmim,n,lvmim,n,zero,emmatim,lda)   !Eim-=(Vim-)^T*sigma_im-
      emmat = emmatre + emmatim
!      
      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvmre,n,zero,smatre,lda)    !Sre=(Vre-)^T*tau_re-
      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvmim,n,zero,smatim,lda)    !Sim=(Vim-)^T*tau_im-
      smat = smatre + smatim
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
      do i_eig = 1, n_max
!
!       being the evec defined as Xeig   = (Yr Zr Yi -Zi)
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
                               dot_product(rmim(:,i_eig),rmim(:,i_eig)))/(sqrt(two)*sqrt2n)
        r_norm(2,i_eig) = maxval(abs(rpre(:,i_eig) + rmre(:,i_eig) + & 
                                     rpim(:,i_eig) + rmim(:,i_eig)))/sqrt(two)
      end do
      !tnwrite(6,*) 'RESIDUAL OLSEN'
      !tndo i_eig = 1, n_max
      !tn  write(6,'(2a20)') 'R+re'//'R+im', 'R-re'//'R-im' 
      !tn  do i = 1, n
      !tn    write(6,'(2f20.10)') rpre(i,i_eig)+rmre(i,i_eig), rpim(i,i_eig)+rmim(i,i_eig)
      !tn  end do
      !tnend do
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
        call ortho_vs_x_complex(.false.,n,ldu,n_act,vpre,vpim,vpre(1,i_beg),vpim(1,i_beg),xx,xx)
        call ortho_vs_x_complex(.false.,n,ldu,n_act,vmre,vmim,vmre(1,i_beg),vmim(1,i_beg),xx,xx)
!        
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
!
!       put current eigenvectors into the first position of the 
!       expansion space
!
        do i_eig = 1, n_max, 2
!         real 
          vpre(:,i_eig)   = (evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
          vmre(:,i_eig)   = (evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
          vpre(:,i_eig+1) = (evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
          vmre(:,i_eig+1) = (evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!         imaginary 
          vpim(:,i_eig)   = ( evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
          vmim(:,i_eig)   = ( evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
          vpim(:,i_eig+1) = (-evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
          vmim(:,i_eig+1) = (-evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
        end do
!
        call ortho_cd_complex(.false.,n,n_max,vpre,vpim,xx,ok)
        call ortho_cd_complex(.false.,n,n_max,vmre,vmim,xx,ok)
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
               bpre,bmre,bpim,bmim)
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
  subroutine caslr_complex_eff_driver(verbose,n,n2,n_targ,n_max,max_iter,tol,max_dav, &
                                         apbmul,ambmul,spdmul,smdmul,lrprec,eig,evecre,evecim,ok)
!
! The following driver is based on the recent SMO-GD algorithm, 
! which efficiently solves the linear response equations, 
! usually defined as:
!
!   / A   B \ / Y \     /  S   D \ / Y \ 
!   |       | |   | = w |        | |   |
!   \ B*  A*/ \ Z*/     \ -D* -S*/ \ Z*/
!
! where A, B, S are hermitian matrices and D is antisymmetric. 
! 
! Although, this routine solves the "swapped" problem, which reads: 
!
!   /  S   D \ / Y \   1 / A   B \ / Y \ 
!   |        | |   | = - |       | |   |
!   \ -D* -S*/ \ Z*/   w \ B*  A*/ \ Z*/ 
!
! For further information read: J. Chem. Theory Comput. 2023, 19, 24, 9025–9031 
!
! Here, we want to solve the complex linear response equations
! using a real algebra formalism.
! The main idea is to apply a sym-antisym separation to each object
! (matrices and vectors) in order to remove the complex algebra, 
! and therefore obtaining a 4nx4n matrix-problem: 
! 
! /  Sr  Dr -Si -Di \/  Yr \       /  Ar  Br -Ai -Bi \ /  Yr \
! | -Dr -Sr -Di -Si ||  Zr | = 1/w |  Br  Ar  Bi  Ai | |  Zr |
! |  Si  Di  Sr  Dr ||  Yi |       |  Ai  Bi  Ar  Br | |  Yi |
! \  Di  Si -Dr -Sr /\ -Zi /       \ -Bi -Ai  Br  Ar / \ -Zi /
!
! As a consequence, the subspace vectors are so defined: 
! 
! b(+) \in R^4n = (b_r+  b_r+ b_i+ -b_i+),
! b(-) \in R^4n = (b_r- -b_r- b_i-  b_i-).
!
! This formalism ensures the orthogonality of the expansion vectors.
!
! Solving the augmented form of the LR-eqs leads to degenerate coupled-eigenvalues, 
! which belong to different eigenvector: 
!
! X1 = ( Yr Zr Yi -Zi), 
! X2 = (-Yi Zi Yr  Zr).
!
! As was already done in the "REAL" solver (caslr_driver) we only manage to deal with
! a reduced form of the expansion vectors; 
! hence, we use only b_r+/b_i+, b_r-/b_i- to build the expansion subspace.
!
! Thanks to the SMO-GD algorithm we build the expansion vectors so that they result
! B-orthogonal (B being the metric).
! So the actual problem we want to solve is a standard eigenvalue problem:
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
!   The first entry of the matrix-vector product subroutines must be:
!   0: real part 
!   1: imaginary part 
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
!   evecre/   double precision real array of size (n2,n_max).
!   evecim:   in input, a guess for the eigenvectors.
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
    real(dp)              :: sqrtn, sqrt2n, tol_rms, tol_max, fac
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
!   eigenvectors of the reduced problem and components of the ritz vectors:
!
!   real
    real(dp), allocatable :: up(:,:), um(:,:), eigpre(:,:), eigmre(:,:), bpre(:,:), bmre(:,:)
!   imaginary
    real(dp), allocatable :: eigpim(:,:), eigmim(:,:), bpim(:,:), bmim(:,:)
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
              bmre(n,n_max), bmim(n,n_max), stat=istat) 
    call check_mem(istat)
!
!   set the tolerances and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    sqrt2n  = sqrt(real(n2,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
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
      vpre(:,i_eig+1) = (evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vmre(:,i_eig+1) = (evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
!     imaginary 
      vpim(:,i_eig)   = ( evecim(1:n,i_eig) - evecim(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig)   = ( evecim(1:n,i_eig) + evecim(n+1:n2,i_eig))/2.0_dp
      vpim(:,i_eig+1) = (-evecre(1:n,i_eig) + evecre(n+1:n2,i_eig))/2.0_dp
      vmim(:,i_eig+1) = (-evecre(1:n,i_eig) - evecre(n+1:n2,i_eig))/2.0_dp
    end do
!
!   orthogonalize the expansion space to the metric.
!
    call apbmul(0,n,n_max,vpre,vpim,lvpre)
    call apbmul(1,n,n_max,vpre,vpim,lvpim)
!    
    call b_ortho_complex(n,n_max,vpre,vpim,lvpre,lvpim)
    
    call ambmul(0,n,n_max,vmre,vmim,lvmre)
    call ambmul(1,n,n_max,vmre,vmim,lvmim)
!    
    call b_ortho_complex(n,n_max,vmre,vmim,lvmre,lvmim)
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
      call spdmul(0,n,n_act,vpre(1,i_beg),vpim(1,i_beg),bvmre(1,i_beg))   !tau_re-
      call spdmul(1,n,n_act,vpre(1,i_beg),vpim(1,i_beg),bvmim(1,i_beg))   !tau_im-
      call smdmul(0,n,n_act,vmre(1,i_beg),vmim(1,i_beg),bvpre(1,i_beg))   !tau_re+
      call smdmul(1,n,n_act,vmre(1,i_beg),vmim(1,i_beg),bvpim(1,i_beg))   !tau_im+
!
      call get_time(t2)
      t_mv = t_mv + t2 - t1
!
!     update the reduced matrix 
!
      call dgemm('t','n',ldu,ldu,n,one,vmre,n,bvmre,n,zero,smatre,lda)    !Sre=(Vre-)^T*tau_re-
      call dgemm('t','n',ldu,ldu,n,one,vmim,n,bvmim,n,zero,smatim,lda)    !Sim=(Vim-)^T*tau_im-
      smat = smatre + smatim
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
!
      do i_eig = 1, n_max
        um(1:ldu,i_eig) = um(1:ldu,i_eig)/eig(i_eig)
      end do
!
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     assemble the symmetric and antysimmetric combinations (Y+Z) and (Y-Z)
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
!     assemble the current approximation to the eigenvectors
!
      do i_eig = 1, n_max
!
!       being the evec defined as Xeig   = (Yr Zr Yi -Zi)
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
        r_norm(1,i_eig) = sqrt(dot_product(rpre(:,i_eig),rpre(:,i_eig)) + &
                               dot_product(rmre(:,i_eig),rmre(:,i_eig)) + &
                               dot_product(rpim(:,i_eig),rpim(:,i_eig)) + &
                               dot_product(rmim(:,i_eig),rmim(:,i_eig)))/(eig(i_eig)*two*sqrt2n)
        r_norm(2,i_eig) = maxval(abs(rpre(:,i_eig) + rmre(:,i_eig) + & 
                                     rpim(:,i_eig) + rmim(:,i_eig)))/(eig(i_eig)*two)
      end do
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
!       check coupled-eigenvalues convergece: if one eigval from the same couple
!       is false, both need to be put false (even if the other was .true.)
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
        call b_ortho_vs_x_complex(n,ldu,n_act,vpre,vpim,lvpre,lvpim,vpre(1,i_beg),vpim(1,i_beg))
!        
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!
        call get_time(t1)
!        
        call apbmul(0,n,n_act,vpre(1,i_beg),vpim(1,i_beg),lvpre(1,i_beg))
        call apbmul(1,n,n_act,vpre(1,i_beg),vpim(1,i_beg),lvpim(1,i_beg))
!        
        call get_time(t2)
        t_mv = t_mv + t2 - t1
!        
        call get_time(t1) 
!        
        call b_ortho_complex(n,n_act,vpre(1,i_beg),vpim(1,i_beg),lvpre(1,i_beg),lvpim(1,i_beg))
!        
        call get_time(t2)
        t_diag = t_diag + t2 - t1
!        
        call get_time(t1) 
!       
        call b_ortho_vs_x_complex(n,ldu,n_act,vmre,vmim,lvmre,lvmim,vmre(1,i_beg),vmim(1,i_beg))
!        
        call get_time(t2)
        t_ortho = t_ortho + t2 - t1
!        
        call get_time(t1) 
!        
        call ambmul(0,n,n_act,vmre(1,i_beg),vmim(1,i_beg),lvmre(1,i_beg))
        call ambmul(1,n,n_act,vmre(1,i_beg),vmim(1,i_beg),lvmim(1,i_beg))
!        
        call get_time(t2)
        t_mv = t_mv + t2 - t1
!        
        call get_time(t1) 
!        
        call b_ortho_complex(n,n_act,vmre(1,i_beg),vmim(1,i_beg),lvmre(1,i_beg),lvmim(1,i_beg))
!        
        call get_time(t2)
        t_diag = t_diag + t2 - t1
!        
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
        call b_ortho_complex(n,n_max,vpre,vpim,lvpre,lvpim)
!        
        call ambmul(0,n,n_max,vmre,vmim,lvmre)
        call ambmul(1,n,n_max,vmre,vmim,lvmim)
!       
        call b_ortho_complex(n,n_max,vmre,vmim,lvmre,lvmim)
!
        bvpre   = zero
        bvmre   = zero
        bvpim   = zero
        bvmim   = zero
        s_red   = zero
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
  subroutine ortho_complex(do_other,n,m,ur,ui,w)
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
    real(dp),  dimension(n,m),   intent(inout) :: ur, ui, w
!
!   local scratch
!   =============
!
    real(dp), allocatable :: vr(:,:), vi(:,:)
!
!   external functions:
!   ===================
!
    external dgeqrf, dtrsm
!
    allocate (vr(n,m),vi(n,m))
    vr = ur
    call dgeqrf(n,m,ur,n,tau,work,lwork,info)
    call dtrsm('r','u','n','n',n,m,one,ur,n,vr,n)
!    
    vi = ui
    call dgeqrf(n,m,ui,n,tau,work,lwork,info)
    call dtrsm('r','u','n','n',n,m,one,ui,n,vi,n)
!
    if (do_other) call dtrsm('r','u','n','n',n,m,one,ur,n,w,n)
!
    ur = vr
    ui = vi
!
    deallocate (vr,vi)
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
  subroutine b_ortho_complex(n,m,ur,ui,bur,bui)
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
    real(dp),  dimension(n,m),   intent(inout) :: ur, ui, bur, bui
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
    call dgemm('t','n',m,m,n,one,ur,n,bur,n,zero,metric,m)
    call dgemm('t','n',m,m,n,one,ui,n,bui,n,one,metric,m)
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
!     TODO
      call dgemm('n','n',n,m,m,one,ur,n,metric,m,zero,temp,n)
      ur = temp
      call dgemm('n','n',n,m,m,one,ui,n,metric,m,zero,temp,n)
      ui = temp
      call dgemm('n','n',n,m,m,one,bur,n,metric,m,zero,temp,n)
      bur = temp
      call dgemm('n','n',n,m,m,one,bui,n,metric,m,zero,temp,n)
      bui = temp
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
      call dtrsm('r','l','t','n',n,m,one,metric,m,ur,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,bur,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,ui,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,bui,n)
    end if
!
    deallocate (metric)
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
  subroutine ortho_cd_complex(do_other,n,m,ur,ui,w,ok)
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
    real(dp),  dimension(n,m), intent(inout) :: ur,ui, w
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
      call dgemm('t','n',m,m,n,one,ur,n,ur,n,zero,metric,m)
      call dgemm('t','n',m,m,n,one,ui,n,ui,n,one,metric,m)
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
        unorm      = dnrm2(n*m,ur,1) + dnrm2(n*m,ui,1)
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
      call dtrsm('r','l','t','n',n,m,one,metric,m,ur,n)
      call dtrsm('r','l','t','n',n,m,one,metric,m,ui,n)
      if (do_other) call dtrsm('r','l','t','n',n,m,one,metric,m,w,n)
!
!     check that the vectors in v are really orthogonal
!
      call dgemm('t','n',m,m,n,one,ur,n,ur,n,zero,metric,m)
      call dgemm('t','n',m,m,n,one,ui,n,ui,n,one,metric,m)
      metric_norm = abs(dnrm2(m*m,metric,1) - sqrt(dble(m)))
      macro_done = metric_norm .lt. tol_ortho
    end do
!
    100 format(t3,'ortho_cd_complex failed with the following error:',a)
!
    ok = .true.
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
      ! I put .and. instead of .or.
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
  subroutine ortho_vs_x_complex(do_other,n,m,k,xr,xi,ur,ui,ax,au)
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
    real(dp),  dimension(n,m), intent(in)    :: xr, xi, ax
    real(dp),  dimension(n,k), intent(inout) :: ur, ui, au
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
    if (.not. useqr) call ortho_cd_complex(do_other,n,k,ur,ui,au,ok)
    if (.not. ok .or. useqr) then 
      !call ortho_complex(do_other,n,k,ur,ui,au)
      write(6,*) 'ortho_cd_complex failed! QR not yet implemented'
    end if 
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^t u)
!
      call dgemm('t','n',m,k,n,one,xr,n,ur,n,zero,xu,m)
      call dgemm('t','n',m,k,n,one,xi,n,ui,n,one,xu,m)
      call dgemm('n','n',n,k,m,-one,xr,n,xu,m,one,ur,n)
      call dgemm('n','n',n,k,m,-one,xi,n,xu,m,one,ui,n)
!
      if (do_other) call dgemm('n','n',n,k,m,-one,ax,n,xu,m,one,au,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd_complex(do_other,n,k,ur,ui,au,ok)
      if (.not. ok .or. useqr) then 
        !call ortho_complex(do_other,n,k,ur,ui,au)
        write(6,*) 'ortho_cd_complex failed! QR not yet implemented'
      end if 
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      call dgemm('t','n',m,k,n,one,xr,n,ur,n,zero,xu,m)
      call dgemm('t','n',m,k,n,one,xi,n,ui,n,one,xu,m)
      xu_norm = dnrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of ortho_vs_x_complex'
    end do
!
    deallocate(xu)
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
  subroutine b_ortho_vs_x_complex(n,m,k,xr,xi,bxr,bxi,ur,ui)
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
    real(dp),  dimension(n,m), intent(in)    :: xr, xi, bxr, bxi
    real(dp),  dimension(n,k), intent(inout) :: ur, ui
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
    if (.not. useqr) call ortho_cd_complex(.false.,n,k,ur,ui,xx,ok)
    if (.not. ok .or. useqr) then 
!      call ortho_complex(.false.,n,k,ur,ui,xx)
      write(6,*) 'ortho_cd_complex failed! QR NOT YET IMPLEMENTED.'
    end if 
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (bx^t u)
!
      call dgemm('t','n',m,k,n,one,bxr,n,ur,n,zero,xu,m)
      call dgemm('t','n',m,k,n,one,bxi,n,ui,n,one,xu,m)
      call dgemm('n','n',n,k,m,-one,xr,n,xu,m,one,ur,n)
      call dgemm('n','n',n,k,m,-one,xi,n,xu,m,one,ui,n)
!
!     now, orthonormalize u.
!
      if (.not. useqr) call ortho_cd_complex(.false.,n,k,ur,ui,xx,ok)
      if (.not. ok .or. useqr) then 
!        call ortho_complex(.false.,n,k,ur,ui,xx)
        write(6,*) 'ortho_cd_complex failed! QR NOT YET IMPLEMENTED.'
      end if 
!
!     compute the overlap between the orthonormalized u and x and decide
!     whether the orthogonalization procedure converged.
!
      call dgemm('t','n',m,k,n,one,bxr,n,ur,n,zero,xu,m)
      call dgemm('t','n',m,k,n,one,bxi,n,ui,n,one,xu,m)
      xu_norm = dnrm2(m*k,xu,1)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of b_ortho_vs_x_complex'
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
