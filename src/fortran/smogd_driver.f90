  module mod_smogd_driver
  use dgl_minor_utils
  use dgl_external_interfaces

  implicit none

  contains

  subroutine smogd_driver(n2,n_targ,n_max,apbmul,ambmul, &
                spdmul,smdmul,lrprec,eig,evec,ok, &
                dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                dgl_memory,dgl_memory_unit)
!!# Driver for the efficient solution to the Linear-Response CASSCF problem
!! \begin{equation}
!!    \begin{bmatrix} \begin{pmatrix}
!!        A & B \\
!!        B & A
!!    \end{pmatrix}
!!    -
!!    \omega
!!    \begin{pmatrix}
!!        S & D \\
!!        -D & -S
!!    \end{pmatrix} \end{bmatrix}
!!    \begin{pmatrix}
!!        Y \\ 
!!        Z
!!    \end{pmatrix}
!!    =
!!    \begin{pmatrix}
!!        0 \\
!!        0
!!    \end{pmatrix},
!!    \label{eq:respeq}
!! \end{equation}
!!
!! Where A, B, S are symmetric matrices and D is antysimmetric.
!!
!! If \(\begin{bmatrix} w, \begin{pmatrix} Y \\ Z \end{pmatrix} \end{bmatrix}\) are a solution,
!! then \(\begin{bmatrix} -w, \begin{pmatrix} Z \\ Y \end{pmatrix} \end{bmatrix}\) is also a solution.
!!
!! Following J. Chem. Phys., 118, 522 (2003), we enforce this property in 
!! the iterative procedure by expanding the eigenvector as
!!
!! \begin{equation}
!! \begin{pmatrix} Y \\ Z \end{pmatrix} =
!! \begin{pmatrix} b^+ \\ b^+ \end{pmatrix} +
!! \begin{pmatrix} b^- \\ -b^- \end{pmatrix}
!! \end{equation}
!!
!! This routine performs the Swapped Metric-Orthogonal -- Generalized Davidsion,
!! therefore solves the associate problem:
!!
!!\begin{equation}
!!    \begin{bmatrix}
!!    \begin{pmatrix}
!!        S & D \\
!!        -D & -S
!!    \end{pmatrix}
!!    -
!!    \frac{1}{\omega}
!!    \begin{pmatrix}
!!        A & B \\
!!        B & A
!!    \end{pmatrix}
!!    \end{bmatrix}
!!    \begin{pmatrix}
!!        Y \\ 
!!        Z
!!    \end{pmatrix}
!!    =
!!    \begin{pmatrix}
!!        0 \\
!!        0
!!    \end{pmatrix},
!!    \label{eq:respeq_smogd}
!!\end{equation}
!!
!! using the casida matrix, which is symmetric and positive definite, as 
!! the metric. This allows us to use expansion vectors that are orthogonal
!! with respect to the dot product defined by the metric, which in turn
!! results in a Rayleigh-Ritz procedure that requires the solution of a 
!! symmetric standard eigenvalue problem
!!
!!\begin{equation}
!!    \begin{pmatrix}
!!        0 & s^T \\
!!        s & 0
!!    \end{pmatrix}
!!    \begin{pmatrix}
!!        u^+ \\ 
!!        u^-
!!    \end{pmatrix}
!!    =
!!    \frac{1}{\omega}
!!    \begin{pmatrix}
!!        u^+ \\ 
!!        u^-
!!    \end{pmatrix},
!!    \label{eq:krylov}
!!\end{equation}
!!
!! which can be reduced to a half-sized eigenvalue problem
!!
!! \(s^T s u^+ = \left(\frac{1}{\omega}\right)^2 u^+ \\
!! u^- = \frac{1}{\omega} Su^+\)
!!
!! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
    implicit none
    integer,                       intent(in)    :: n2
!! Totals size of the generalized eigenvalue
    integer,                       intent(in)    :: n_targ
!! Number of required eigenpairs.
    integer,                       intent(in)    :: n_max
!! Maximum size of the search space. Should be >= n_targ.
    real(dp), dimension(n_max),    intent(inout) :: eig
!! Computed eigenvalues
    real(dp), dimension(n2,n_max), intent(inout) :: evec
!! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
    logical,                       intent(inout) :: ok
!! True if davidson converged
    procedure(smogd_matvec) :: apbmul
!! External subroutine that performs the matrix-vector multiplication with A+B
    procedure(smogd_matvec) :: ambmul
!! External subroutine that performs the matrix-vector multiplication with A-B
    procedure(smogd_matvec) :: spdmul
!! External subroutine that performs the matrix-vector multiplication with S+D
    procedure(smogd_matvec) :: smdmul
!! External subroutine that performs the matrix-vector multiplication with S-D
    procedure(smogd_precnd) :: lrprec
!! External subroutine that applies a preconditioner to both plus and minus vectors
    logical,  optional,            intent(in)    :: dgl_verbose
!! Verbose mode. Default = .false.
    integer,  optional,            intent(in)    :: dgl_dav_iter
!! Maximum number of iterations before Davidson restart. Default = \(25\)
    integer,  optional,            intent(in)    :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
    integer,  optional,            intent(in)    :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
    character(len=2),  optional,   intent(in)    :: dgl_memory_unit
!! Unit of memory. Default = MBs
    real(dp), optional,            intent(in)    :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
!
!   local variables:
!   ================
    logical  :: bool
    integer  :: max_iter, dav_iter, memory
    real(dp) :: tol
    character(len=2) :: memory_unit
!
!   dimension of the halved problem, the one we are actually solving
!
    integer               :: n
!
!   actual expansion space size and total dimension
!
    integer               :: lda, lda2
!
!   number of active vectors at a given iteration, and indices to access them
!
    integer               :: n_act, ind, i_beg
!
!   current size and total dimension of the expansion space
!
    integer               :: ld_current
!
!   number of large arrays that will be allocated
!
    integer               :: n_arrs
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
    real(dp), allocatable :: s_copy(:,:), s_red_2(:,:), e_red(:)
    real(dp), allocatable :: s_red(:,:)
!
!   Scratch vector to avoid recomputation of reduced matrix
!
    real(dp), allocatable :: scratch(:,:)
!
!   ================
!   START EXECUTION
!   ================
!
!   Stupidity checks
!
    if(n_targ.gt.n_max) call dgl_error(&
    "Number of eigenvalues request is larger that size of arrays passed")
    if(mod(n2,2).ne.0) call dgl_error(&
    "Size of the total problem is not even, something is really wrong with your input")
!
!   Parse optional arguments
!
    bool = .false.     ; if(present(dgl_verbose)) bool = dgl_verbose
    max_iter = 100     ; if(present(dgl_max_iter)) max_iter = dgl_max_iter
    dav_iter = 25      ; if(present(dgl_dav_iter)) dav_iter = dgl_dav_iter
    tol = 1.e-7_dp     ; if(present(dgl_tol)) tol = dgl_tol
    memory = 80        ; if(present(dgl_memory)) memory = dgl_memory
    memory_unit = "MB" ; if(present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
!
!   compute the actual size of the expansion space, checking that
!   the input makes sense.
!   no expansion space smaller than dav_iter = 10 is deemed acceptable.
!
    n    = n2 / 2
    lda  = dav_iter*n_max
    lda2 = 2 * lda
    if (lda .ge. n) then
      if (bool) call dgl_warning("Expansion space is larger than the dimension of the problem. " //&
                       "Reducing size to avoid Rouché-Capelli failure" )
      lda = n - 1
    endif
!
!   compute the number of large vectors that will be allocated 
!   to later exstimate required memory in dgl_init 
!
    n_arrs = lda*6 + n_max*7
    call dgl_init(n,n_arrs,memory,memory_unit,bool)
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max)
    call mallocate (lwork,work)
    call mallocate (n_max,tau)
!
!   allocate memory for the expansion space, the corresponding 
!   matrix-multiplied vectors and the residual:
!
    call mallocate(n,lda,vp)
    call mallocate(n,lda,vm)
    call mallocate(n,lda,lvp)
    call mallocate(n,lda,lvm)
    call mallocate(n,lda,bvp)
    call mallocate(n,lda,bvm)
    call mallocate(n,n_max,rp)
    call mallocate(n,n_max,rm)
    call mallocate(n,n_max,rr)
!
!   allocate memory for convergence check
!
    call mallocate(n_max,done)
    call mallocate(2,n_max,r_norm)
!
!   allocate memory for the reduced matrix and its eigenvalues:
!
    call mallocate(lda,lda,s_copy)
    call mallocate(lda,lda,s_red_2)
    call mallocate(lda,lda,s_red)
    call mallocate(lda2,e_red)
!
!   allocate memory for the plus and minus eigenvector components:
!
    call mallocate(lda,n_max,up)
    call mallocate(lda,n_max,um)
    call mallocate(n,n_max,eigp)
    call mallocate(n,n_max,eigm)
    call mallocate(n,n_max,bp)
    call mallocate(n,n_max,bm)
!
    call mallocate(n_max,lda,scratch)
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
!   initialize the counters
!
    n_act = n_max
    ind   = 1
    i_beg = 1
    ld_current = 0
!
!   main loop:
!
    1030 format(t5,'SMO-GD iterations (tol=',d10.2,'):',/, &
                t5,'------------------------------------------------------------------',/, &
                t7,'  iter  root              eigenvalue','         rms         max ok',/, &
                t5,'------------------------------------------------------------------')
    1040 format(t9,i4,2x,i4,f24.12,2d12.4,l3)
!
    if (verbose) write(6,1030) tol
!
    do it = 1, max_iter
!
!     update the size of the expansion space.
!
      ld_current = ld_current + n_act
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
      call dgemm('t','n',ld_current,n_act,n,one,vm,n,bvm(:,i_beg),n,zero,s_red(1,i_beg),lda)
      if (it.gt.1) then
        call dgemm('t','n',n_act,i_beg-1,n,one,vm(:,i_beg),n,bvm,n,zero,scratch,n_max)
        s_red(i_beg:ld_current,1:i_beg-1) = scratch(:n_act,:i_beg-1)
      endif
!
!     save s, and assemble s^t s:
!
      s_copy  = s_red
      call dgemm('t','n',ld_current,ld_current,ld_current,one,s_copy,lda,s_copy,lda,zero,s_red_2,lda)
!
!     diagonalize s^t s
!
      call get_time(t1)
      call dsyev('v','u',ld_current,s_red_2,lda,e_red,work,lwork,info)
      call get_time(t2)
      t_diag = t_diag + t2 - t1
!
!     extract the eigenvalues and compute the ritz approximation to the
!     eigenvectors 
!
      do i_eig = 1, n_max
        eig(i_eig)      = sqrt(e_red(ld_current - i_eig + 1))
        up(1:ld_current,i_eig) = s_red_2(1:ld_current,ld_current - i_eig + 1)
      end do
!
!     compute the u_- eigenvectors:
!
      call dgemm('n','n',ld_current,n_max,ld_current,one,s_copy,lda,up,lda,zero,um,lda)
      do i_eig = 1, n_max
        um(1:ld_current,i_eig) = um(1:ld_current,i_eig)/eig(i_eig)
      end do
!
!     asemble the symmetric and antysimmetric combinations (Y+Z) and (Y-Z)
!
      call dgemm('n','n',n,n_max,ld_current,one,vp,n,up,lda,zero,eigp,n)
      call dgemm('n','n',n,n_max,ld_current,one,vm,n,um,lda,zero,eigm,n)
!
!     assemble the current approximation to the eigenvectors
!
      evec(1:n,:)    = eigp + eigm
      evec(n+1:n2,:) = eigp - eigm
!
!     compute the residuals, and their rms and sup norms:
!
      call dgemm('n','n',n,n_max,ld_current,one,bvp,n,um,lda,zero,rp,n)
      call dgemm('n','n',n,n_max,ld_current,one,bvm,n,up,lda,zero,rm,n)
      call dgemm('n','n',n,n_max,ld_current,one,lvp,n,up,lda,zero,bp,n)
      call dgemm('n','n',n,n_max,ld_current,one,lvm,n,um,lda,zero,bm,n)
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
      if (ld_current + n_act .le. lda) then
!
        i_beg = i_beg + n_act
!        
      else
!
        if (verbose) write(6,'(t7,a)') 'Restarting davidson.'
!        
!       put current eigenvectors into the first position of the 
!       expansion space
!
        vp(:,:n_max) = eigp
        vm(:,:n_max) = eigm
!
        lvp(:,:n_max) = bp
        lvm(:,:n_max) = bm
        call b_ortho(n,n_max,vp,lvp)
        call b_ortho(n,n_max,vm,lvm)
!
        call dgemm('n','n',n,n_max,ld_current,one,bvp,n,um,lda,zero,bp,n)
        call dgemm('n','n',n,n_max,ld_current,one,bvm,n,up,lda,zero,bm,n)
        bvp(:,:n_max) = bp
        bvm(:,:n_max) = bm
!
        s_red = zero
        do i_eig = 1, n_max
          s_red(i_eig,i_eig) = eig(i_eig)
        enddo
!
!       initialize indexes back to their starting values 
!
        ld_current = n_max
        i_beg = n_max + 1
!
      end if
!
!       compute the preconditioned residuals using davidson's procedure
!       note that this is done with a user-supplied subroutine, that can
!       be generalized to experiment with fancy preconditioners that may
!       be more effective than the diagonal one, as in the original 
!       algorithm.
!
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
!     orthogonalize the new vectors to the existing ones and then
!     orthonormalize them.
!
      call get_time(t1)
      call b_ortho_vs_x(n,ld_current,n_act,vp,lvp,vp(1,i_beg))
      call apbmul(n,n_act,vp(1,i_beg),lvp(1,i_beg))
      call b_ortho(n,n_act,vp(1,i_beg),lvp(1,i_beg))
      call b_ortho_vs_x(n,ld_current,n_act,vm,lvm,vm(1,i_beg))
      call ambmul(n,n_act,vm(1,i_beg),lvm(1,i_beg))
      call b_ortho(n,n_act,vm(1,i_beg),lvm(1,i_beg))
      call get_time(t2)
      t_ortho = t_ortho + t2 - t1
!      
      if (verbose) write(6,1050) n_targ, n_act, n_frozen
!      
    end do
!
    call get_time(t1)
    t_tot = t1 - t_tot
!    
    if (verbose) write(6,1000) t_mv, t_diag, t_ortho, t_tot
!    
    call mfree(work)
    call mfree(tau)
    call mfree(vp)
    call mfree(vm)
    call mfree(lvp)
    call mfree(lvm)
    call mfree(bvp)
    call mfree(bvm)
    call mfree(rp)
    call mfree(rm)
    call mfree(rr)
    call mfree(r_norm)
    call mfree(done)
    call mfree(s_copy)
    call mfree(s_red_2)
    call mfree(s_red)
    call mfree(e_red)
    call mfree(up)
    call mfree(um)
    call mfree(eigp)
    call mfree(eigm)
    call mfree(bp)
    call mfree(bm)
    call mfree(scratch)
!
    call dgl_check_memleak()
!
    1000 format(t3,'timings for caslr_eff (cpu/wall):   ',/, &
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
  end subroutine smogd_driver

end module mod_smogd_driver