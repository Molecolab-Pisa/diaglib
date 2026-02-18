module mod_davidson_nosym_driver
  use dgl_minor_utils
  use dgl_external_interfaces

implicit none

contains

subroutine davidson_nosym_driver(n,n_targ,n_max,matvec_r,matvec_l,precnd,side, &
              eig,evec_r,evec_l,ok,&
              dgl_verbose,dgl_tol,dgl_max_iter,dgl_dav_iter,&
              dgl_shift,dgl_memory,dgl_memory_unit)
!! ### Driver for Davidson-Liu non-symmetric diagonalization
!! Can solve solve only standard eigenvalue problems. Can eveluate both Left and Right eigenvectors.
!! To pass any optional argument, you need to use [[dgl_drivers_interfaces]], a module included in this library.
!!
!! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
  implicit none
    integer,                      intent(in)    :: n
!! Size of the matrix to be diagonalized
    integer,                      intent(in)    :: n_targ
!! Number of required eigenpairs.
    integer,                      intent(in)    :: n_max
!! Maximum size of the search space. Should be >= n_targ
    integer,                      intent(in)    :: side
!! Integer to decide which eigenvectors to compute and whether to compute
!! them togheter or separately
    real(dp), dimension(n_max),   intent(inout) :: eig
!! Computed eigenvalues
    real(dp), dimension(n,n_max), intent(inout) :: evec_l
!! Computed Left eigenvectors. In input, it should contain their guess
    real(dp), dimension(n,n_max), intent(inout) :: evec_r
!! Computed Right eigenvectors. In input, it should contain their guess
    logical,                      intent(inout) :: ok
!! True if davidson converged
    procedure(matvec_) :: matvec_r
!! External subroutine that performs the matrix-vector multiplication for
!! right eigenvectors
    procedure(matvec_) :: matvec_l
!! External subroutine that performs the matrix-vector multiplication for
!! left eigenvectors
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
!
!   local variables:
!   ================
    logical  :: verbose
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
!   restarting variables
!
    logical               :: restart
!
!   iterators and utilities
!
    integer               :: it, i_eig
    real(dp)              :: sqrtn, tol_im
    real(dp)              :: xx(1), yy   
    integer               :: j, k
!
!   arrays to control convergence
!
    logical, allocatable  :: done(:)
!
!   expansion spaces, residuals and their norms
!
    real(dp), allocatable :: space_r(:,:), space_l(:,:), aspace_r(:,:), aspace_l(:,:)
    real(dp), allocatable :: residuals_r(:,:), residuals_l(:,:)
    real(dp), allocatable :: r_norm_r(:,:), r_norm_l(:,:)
!
!   subspace matrix, eigenvalues and real and imaginary parts of the eigenvalues
!
    real(dp), allocatable :: a_red(:,:), e_red_re(:), e_red_im(:)
    real(dp), allocatable :: evec_red_r(:,:), evec_red_l(:,:)
    real(dp), allocatable :: copy_r(:,:), copy_l(:,:), copy_eig(:)
!
!   variables for left, right, or both eigenvectors
!
    logical               :: left, right, consecutive, do_davidson
    real(dp)              :: eig_r(n_max)
! 
!   variables for the sorting eigenpairs, since lapack does not.
!
    integer               :: max_idx(1)
    logical               :: found_im, found_er, double_r, double_l, ortho_ok
    logical, allocatable  :: mask_overlap(:)
!
    real(dp),allocatable  :: overlap(:,:), perm_mat(:,:), evec_temp(:,:), eig_temp(:), overlap_diff(:)
    real(dp)              :: overlap_idx_r(n_max,2), overlap_val_r(n_max,2), overlap_self_l(n_max), &
                             overlap_idx_l(n_max,2), overlap_val_l(n_max,2), overlap_self_r(n_max)
    real(dp),allocatable  :: perm_temp(:,:)
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
! Parse optional arguments
!
    verbose = .false.   ; if(present(dgl_verbose)) verbose = dgl_verbose
    max_iter = 100      ; if(present(dgl_max_iter)) max_iter = dgl_max_iter
    dav_iter = 25       ; if(present(dgl_dav_iter)) dav_iter = dgl_dav_iter
    tol = 1.e-7_dp      ; if(present(dgl_tol)) tol = dgl_tol
    shift = 0.e0_dp     ; if(present(dgl_shift)) shift = dgl_shift
    memory = 80         ; if(present(dgl_memory)) memory = dgl_memory
    memory_unit = "MB"  ; if(present(dgl_memory_unit)) memory_unit = dgl_memory_unit
!
!   set some quantities
!
    ok          = .false.
    right       = .false.
    left        = .false.
    consecutive = .false.
    do_davidson = .true.
!
!   extract from the input which eigenvectors shall be computed in which way
!     1 = only right eigenpairs
!     2 = only left eigenpairs
!     3 = both eigenpairs in simultaneous manner
!     4 = both eigenpairs in consecutive manner, start with right
!
    if (side .eq. 1) then 
      right = .true.
    else if (side .eq. 2) then 
      left  = .true.
    else if (side .eq. 3) then
!
!     the simultaneous diagonalization driver is less efficient than a consecutive
!     run for the left eigenvectors. 
!     switch it off manually, but leave it as an advanced debug feature.
!
      right = .true.
      consecutive = .true.
    else if (side .eq. 4) then 
      right = .true.
      consecutive = .true.
    else
      call dgl_error(" Choice for side is not correct. " // achar(10) // &
      "   Has to be 1(right), 2(left), 3(both simultaneuos) or 4(both consequentially).")
    end if
!
!   computing actual size of the expansion space, checking that
!   the input makes sense.
!
    lda = dav_iter*n_max
!
!   compute the number of large vectors that will be allocated 
!   to later exstimate required memory in dgl_init 
!
    n_arrs = lda*4 + n_max*2
    call dgl_init(n,n_arrs,memory,memory_unit,verbose)
!
!   start by allocating memory for the various lapack routines
!
    lwork = get_mem_lapack(n,n_max) 
    call mallocate(lwork,work)
    call mallocate(lda,tau)
!
!   allocate memory for for expansion space, the corresponding
!   matrix-multiplied vectors and the residuals
!
    call mallocate(n,lda,space_r)
    call mallocate(n,lda,space_l)
    call mallocate(n,lda,aspace_l)
    call mallocate(n,lda,aspace_r)
    call mallocate(n,n_max,residuals_l)
    call mallocate(n,n_max,residuals_r)
!
!   allocate memory for convergency check
!
    call mallocate(n_max,done)
    call mallocate(2,n_max,r_norm_l)
    call mallocate(2,n_max,r_norm_r)
!
!   allocate memory for the reduced matrix, its eigenvalues with real &
!   imaginary parts, and its left & right eigenvectors 
!
    call mallocate(lda,lda,a_red)
    call mallocate(2*lda,e_red_re)
    call mallocate(2*lda,e_red_im)
    call mallocate(lda,lda,evec_red_l)
    call mallocate(lda,lda,evec_red_r)
    call mallocate(lda,lda,copy_l)
    call mallocate(lda,lda,copy_r)
    call mallocate(2*lda,copy_eig)
!
!   allocate space for orthogonalization routines
!   and mask array for sorting routine
!
    call mallocate(2*n_max,2*n_max,overlap)
    call mallocate(n_max,overlap_diff)
    call mallocate(2*n_max,2*n_max,perm_mat)
    call mallocate(2*n_max,2*n_max,perm_temp)
    call mallocate(lda,n_max,evec_temp)
    call mallocate(lda,eig_temp)
    call mallocate(2*n_max,mask_overlap)
!
!   set the tolerance and compute a useful constant to compute rms norms:
!
    sqrtn   = sqrt(real(n,dp))
    tol_rms = tol
    tol_max = 10.0_dp * tol
    tol_im  = 1.d-12
!
!   check weather we have a guess for the eigenvectors in evec, and
!   weather it is orthonormal.
!   if evec is zero, create a random guess
!
    if (right) call check_guess(n,n_max,evec_r)
    if (left)  call check_guess(n,n_max,evec_l)
!
    do while (do_davidson)
! 
!     clean out various quantities
!
      t_tot       = zero
      t_diag      = zero
      t_ortho     = zero
      t_mv        = zero
      space_r     = zero
      space_l     = zero
      aspace_r    = zero
      aspace_l    = zero
      a_red       = zero
      e_red_re    = zero
      e_red_im    = zero
      copy_r      = zero
      copy_l      = zero
      r_norm_r    = zero
      r_norm_l    = zero
      residuals_r = zero
      residuals_l = zero
      done = .false.
!
      call get_time(t_tot1)
!
!     move guess into the expansion spaces
!
      if (right) call dcopy(n*n_max,evec_r,1,space_r,1)
      if (left)  call dcopy(n*n_max,evec_l,1,space_l,1)
!
      n_act = n_max
      i_beg = 1
      ind   = 1
!
!     initialize the counter for the expansion of the subspace
!
      ld_current   = 0
!
!     initialize to false the restart
!
      restart = .false.
!
!     print header
!
      if (verbose) write(6,1030) tol, left, right
! 
!     main loop
!     
      do it = 1, max_iter
!
!       update the size of the expansion space.
!
        ld_current   = ld_current + n_act
!
!       perform this iteration's matrix-vector multiplications for both 
!       right and left expansion spaces
!
        call get_time(t1)
        if (right) call matvec_r(n,n_act,space_r(1,i_beg),aspace_r(1,i_beg))
        if (left)  call matvec_l(n,n_act,space_l(1,i_beg),aspace_l(1,i_beg))
        call get_time(t2)
        t_mv = t_mv + t2 -t1
!
!       get the reduced matrix
!
        if (left .and. right) then 
          call dgemm('t','n',ld_current,ld_current,n,one,space_l,n,aspace_r,n,zero,a_red,lda)
        else if (right) then
          call dgemm('t','n',ld_current,ld_current,n,one,space_r,n,aspace_r,n,zero,a_red,lda)
        else if (left) then
          call dgemm('t','n',ld_current,ld_current,n,one,aspace_l,n,space_l,n,zero,a_red,lda)
        end if
!
!       diagonalize the reduced matrix
!
        call get_time(t1)
        call dgeev('v','v',ld_current,a_red,lda,e_red_re,e_red_im,evec_red_l,lda,evec_red_r,lda,work,lwork,info)
        call get_time(t2)
!
        t_diag = t_diag + t2 - t1
! 
!       check if diagonalization terminated with info = 0
!
        if (info.ne.0) then
          call dgl_error("diagonalization of reduced space failed.")
        end if
!
!       sort lowest eigenpairs in increasing order in range 2*n_max to ensure that all n_max 
!       sought eigenpairs are in the range 2*n_max
!
        if (it.gt.1 .and. .not. restart) then
!
          call sort_eigenpairs(ld_current,e_red_re,e_red_im,evec_red_r,evec_red_l,n_max+n_act,lda,.true.,tol_im)
!
        else if (it.eq.1 .or. restart) then
!
          call sort_eigenpairs(ld_current,e_red_re,e_red_im,evec_red_r,evec_red_l,n_max,lda,.true.,tol_im)
!
        end if
!
!       double check for complex contributions in the n_max sought eigenvalues
!
        found_im = .false.
        do j = 1, n_max
          if (e_red_im(j).gt.tol_im) found_im = .true.
        end do
!
        if (found_im.and.verbose) then
          print *
          print *, "complex contribution in sought eigenvalues"
          print *
        end if
!
!       compute overlap of old and new eigenvectors in the dimension of the old eigenvectors
!       to ensure correct sorting by checking if largest absolute value of column is on the
!       diagonal. if not, use the indices of the largest elements to construct a permutation
!       matrix to resort the eigenpairs according to the overlap.
!  
        if (it.ne.1 .and. .not. restart) then
!
!         compute overlap for the right eigenvectors and extract the index and value of the largest
!         and second largest overlap
!
          call dgemm('t','n',2*n_max,2*n_max,ld_current,one,copy_r,lda,evec_red_r,lda,zero,overlap,2*n_max)
!
          found_er = .false.
          do j = 1, n_max
            mask_overlap = .true.
            max_idx = maxloc(abs(overlap(:,j)))
            overlap_idx_r(j,1)  = real(max_idx(1), kind=dp)
            overlap_self_r(j) = overlap(j,j)
            overlap_val_r(j,1) = overlap(max_idx(1),j)
            mask_overlap(max_idx) = .false.
!           
!           identify if a swapping is necessary
!
            if (max_idx(1).ne.j) then
                found_er = .true.
            end if
!
!           extract index and value of second larges overlap
!
            max_idx = maxloc(abs(overlap(:,j)), mask =mask_overlap)
            overlap_idx_r(j,2)  = real(max_idx(1), kind=dp)
            overlap_val_r(j,2) = overlap(max_idx(1),j)
          end do
!
!         compute overlap for the left eigenvectors
!
          call dgemm('t','n',2*n_max,2*n_max,ld_current,one,copy_l,lda,evec_red_l,lda,zero,overlap,2*n_max)
!
          do j = 1, n_max
            mask_overlap = .true.
            max_idx = maxloc(abs(overlap(:,j)))
            overlap_idx_l(j,1)  = real(max_idx(1), kind=dp)
            overlap_self_l(j) = overlap(j,j)
            overlap_val_l(j,1) = overlap(max_idx(1),j)
            mask_overlap(max_idx) = .false.
!
!           identify if a swapping is necessary
!
            if (max_idx(1).ne.j) then
              found_er = .true.
            end if
!
!           extract index and value of second largest overlap
!
            max_idx = maxloc(abs(overlap(:,j)), mask =mask_overlap)
            overlap_idx_l(j,2)  = real(max_idx(1), kind=dp)
            overlap_val_l(j,2) = overlap(max_idx(1),j)
          end do
!
!         check if no indices were assigned twice as maximum overlap
!
          double_r = .false.
          double_l = .false.
          do j=1, n_max
            do k=1, n_max
              if (k .ne. j .and. abs(overlap_idx_r(j,1) - overlap_idx_r(k,1)) .lt. num_thresh) then 
                double_r = .true.
              end if
            end do
          end do
          do j=1, n_max
            do k=1, n_max
              if (k .ne. j .and. abs(overlap_idx_l(j,1) - overlap_idx_l(k,1)) .lt. num_thresh) then 
                double_l = .true.
              end if
            end do
          end do
!
!         try easy fix, by just taking the permutation idices of the other eigenvector side
!
          if (double_r .and. .not. double_l) then
            overlap_idx_r(:,1) = overlap_idx_l(:,1)
          else if (double_l .and. .not. double_r) then
            overlap_idx_l(:,1) = overlap_idx_r(:,1)
          else if (double_r .and. double_l) then
!
!           check which second largest overlap is larger and take this indice as max_overlap.
!           try for right side only, if no result, dont swap anything and try to continue 
!           without swapping any eigenvectors
!
            do j=1, n_max
              do k=1, n_max
                if (k .ne. j .and. abs(overlap_idx_r(j,1) - overlap_idx_r(k,1)) .lt. num_thresh) then 
                  if (overlap_val_r(j,2) .gt. overlap_val_r(k,2)) then
                    overlap_idx_r(j,1) = overlap_idx_r(j,2)
                  else
                    overlap_idx_r(k,1) = overlap_idx_r(k,2)
                  end if
                end if
              end do
            end do
!
!           check again if they are the same indices in the max_overlap for the right side.
!           if no, then take these indices for the right and left side. if yes, try without 
!           swapping
!
            double_r = .false.
            do j=1, n_max
              do k=1, n_max
                if (k .ne. j .and. abs(overlap_idx_r(j,1) - overlap_idx_r(k,1)) .lt. num_thresh) then 
                  double_r = .true.
                end if
              end do
            end do
            if (double_r) then
              do j=1, n_max
                overlap_idx_r(j,1) = real(j, kind=dp)
                overlap_idx_l(j,1) = real(j, kind=dp)
              end do
            else
              overlap_idx_l(:,1) = overlap_idx_r(:,1)
            end if  
          end if
!
!         check if maximum indices from right and left eigenvectors are the same
!
          overlap_diff = overlap_idx_r(:,1) - overlap_idx_l(:,1)
          if (dnrm2(n_max,overlap_diff,1 ) .gt. num_thresh) then
            if (sum(overlap_val_r(:,1)) .gt. sum(overlap_val_l(:,1))) then
              overlap_idx_l(:,1) = overlap_idx_r(:,1)
            else
              overlap_idx_r(:,1) = overlap_idx_l(:,1)
            end if
          end if
!
          if (found_er) then
!
!             now permute eigenvectors according to the maxiumum overlap.
!             get permutation matrix first
!
              perm_mat = zero 
              do j=1, n_max
                perm_mat(int(overlap_idx_r(j,1)),j) = one
              end do
              perm_temp = transpose(perm_mat)
!
!             now permute left & right eigenvectors and imaginary & real eigenvalues 
!             note: the use of 't' instead of computing the transpose explicitly obtained in
!                   a different result
!
              call dgemm('n','n',ld_current,n_max,2*n_max,one,evec_red_r,lda,perm_mat,2*n_max,zero,evec_temp,lda)
              call dcopy(ld_current*n_max,evec_temp,1,evec_red_r,1)
!
              call dgemm('n','n',ld_current,n_max,2*n_max,one,evec_red_l,lda,perm_mat,2*n_max,zero,evec_temp,lda)
              call dcopy(ld_current*n_max,evec_temp,1,evec_red_l,1)
!
              call dgemv('n',n_max,2*n_max,one,perm_temp,2*n_max,e_red_re,1,zero,eig_temp,1)
              call dcopy(n_max,eig_temp,1,e_red_re,1)
!
              call dgemv('n',n_max,2*n_max,one,perm_temp,2*n_max,e_red_im,1,zero,eig_temp,1)
              call dcopy(n_max,eig_temp,1,e_red_im,1)
!
          end if
        end if
!
!       copy and save the new eigenvectors for the next iteration
!
        restart = .false.
        copy_r   = evec_red_r
        copy_l   = evec_red_l
        copy_eig = e_red_re
!
!       extract the eigenvalues and compute the ritz approximation to the 
!       eigenvectors
!
        eig = e_red_re(1:n_max)
!
        if (right) call dgemm('n','n',n,n_max,ld_current,one,space_r,n,evec_red_r,lda,zero,evec_r,n)
        if (left)  call dgemm('n','n',n,n_max,ld_current,one,space_l,n,evec_red_l,lda,zero,evec_l,n)
!
!       compute the residuals, and their rms and sup norms
!
        if (right) call dgemm('n','n',n,n_max,ld_current,one,aspace_r,n,evec_red_r,lda,zero,residuals_r,n) 
        if (left)  call dgemm('n','n',n,n_max,ld_current,one,aspace_l,n,evec_red_l,lda,zero,residuals_l,n) 
!
        do i_eig = 1, n_targ
!
!         if the eigenvalue is already converged, skip it.
!
          if (done(i_eig)) cycle
!
          if (right) then
            call daxpy(n,-eig(i_eig),evec_r(:,i_eig),1,residuals_r(:,i_eig),1)
            r_norm_r(1,i_eig) = dnrm2(n,residuals_r(:,i_eig),1)/sqrtn
            r_norm_r(2,i_eig) = maxval(abs(residuals_r(:,i_eig)))
          end if
          if (left) then
            call daxpy(n,-eig(i_eig),evec_l(:,i_eig),1,residuals_l(:,i_eig),1)
            r_norm_l(1,i_eig) = dnrm2(n,residuals_l(:,i_eig),1)/sqrtn
            r_norm_l(2,i_eig) = maxval(abs(residuals_l(:,i_eig)))
          end if 
!
        end do
!
!       check convergence. lock the first contiguous converged eigenvalues
!       by setting the logical array "done" to true
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
!       print some information
!
        if (verbose) then
          do i_eig = 1, n_targ
            write(6,1040) it, i_eig, eig(i_eig) - shift, r_norm_l(:,i_eig), r_norm_r(:,i_eig), done(i_eig)
          end do
          write(6,*)
        end if 
!
        if (all(done(1:n_targ))) then
          ok = .true.
          exit
        end if
!       
!       check weather an update is required.
!       if not, perform a davidson restart
!
        if (ld_current + n_act .lt. lda) then 
!
!         compute the preconditioned residuals using davidson's procedure
!         note that this is done with a user-supplied subroutine, that can
!         be generalized to experiment with fancy preconditioners that may
!         be more effective than the diagonal one, as in the original 
!         algorithm.
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
          if (right) call precnd(n,n_act,-eig(ind),residuals_r(1,ind),space_r(1,i_beg))
          if (left)  call precnd(n,n_act,-eig(ind),residuals_l(1,ind),space_l(1,i_beg))
!
!         orthogonalize the new vectors to the existing ones of the respective other
!         space and orthogonalize set of new vectors among each other
!
!         Gram-Schmit orthogonalization of residual to the respective subspace 
!
          call get_time(t1)
          if (right .and. left) then
            call biortho_vs_x(n,ld_current,n_act,space_l,space_r,space_l(1,i_beg),space_r(1,i_beg))
          else if (right) then
            call ortho_vs_x(n,ld_current,n_act,space_r,space_r(1,i_beg),xx,xx)
          else if (left) then
            call ortho_vs_x(n,ld_current,n_act,space_l,space_l(1,i_beg),xx,xx)
          end if
          call get_time(t2)
!
          t_ortho = t_ortho + t2 - t1
!
!         normalize columns 
!
        else 
          if (verbose) write(6,'(t7,a)') 'Restarting Davidson'
          n_act   = n_max
          space_r = zero
          space_l = zero
!
!         put current eigenvectors into the first position of tne
!         expansion space
!
          if (right) call dcopy(n_max*n,evec_r,1,space_r,1)
          if (left)  call dcopy(n_max*n,evec_l,1,space_l,1)
!
          ortho_ok = .false.
          call get_time(t1)
          if (right .and. left) then
            call svd_biortho(n,n_act,space_r,space_l)
          else if (right) then
            call ortho_cd(n,n_act,space_r,yy,ok) 
          else if (left) then
            call ortho_cd(n,n_act,space_l,yy,ok) 
          end if
          call get_time(t2)
          t_ortho = t_ortho + t2 - t1
!
          aspace_r  = zero
          aspace_l  = zero
          a_red     = zero
          e_red_re  = zero
          e_red_im  = zero
!
!         initialize indexes back to their starting values
!
          ld_current   = 0
          i_beg = 1
!
!         counting how many matvec we can skip at the next
!         iteration
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
!     end of davidson, print results
!
      call get_time(t_tot2)
      t_tot = t_tot2 - t_tot1
!
!     if required, print timings
!
      if (verbose) then
        print * 
        write(6,1100) t_mv, t_diag, t_ortho, t_tot
        print * 
        print * 
      end if
!
!     stop after one davidson evaluation or do a second one if side = consecutive 
!
      if (side.eq.1 .or. side.eq.2) do_davidson = .false.
      if (consecutive) then
        if (left) then
          left = .false.
          do_davidson = .false.
!
!         check if energies are same
!
          if (maxval(eig_r(:n_targ) - eig(:n_targ)) .gt. tol) then
            print *, "Debug: eigenvalues in the consecutive computation of", &
                     "right and left eigenpairs do not match. Stopping" 
            stop
          end if
!
        else
          right = .false.
          left  = .true.
          eig_r = eig
!
!         use evec_r as guess for evec_l
!
          call ortho_cd(n,n_max,evec_r,yy,ok) 
          call dcopy(n*n_max,evec_r,1,evec_l,1) 
        end if  
        
      end if
    end do
!
!   final orthogonalization of  evec_r and evec_l, if both were computed.
!
    if (consecutive .or. left .and. right) then
      call get_time(t1)
      call svd_biortho(n,n_max,evec_l,evec_r)
      call get_time(t2)
!
      t_ortho = t_ortho + t2 - t1
    end if
!      
!   deallocate memory
!
    call mfree(work)
    call mfree(tau)
    call mfree(space_r)
    call mfree(space_l)
    call mfree(aspace_l)
    call mfree(aspace_r)
    call mfree(residuals_l)
    call mfree(residuals_r)
    call mfree(r_norm_l)
    call mfree(r_norm_r)
    call mfree(done)
    call mfree(a_red)
    call mfree(e_red_re)
    call mfree(e_red_im)
    call mfree(evec_red_l)
    call mfree(evec_red_r)
    call mfree(copy_l)
    call mfree(copy_r)
    call mfree(copy_eig)
    call mfree(overlap)
    call mfree(overlap_diff)
    call mfree(perm_mat)
    call mfree(perm_temp)
    call mfree(evec_temp)
    call mfree(eig_temp)
    call mfree(mask_overlap)
!
    call dgl_check_memleak()
!
1100 format(t3,'  timings for non-symmetric Davidson (cpu/wall) : ',/, &
            t3,'  matrix-vector multiplications   : ',2f12.4,/, &
            t3,'  diagonalization                 : ',2f12.4,/, &
            t3,'  orthogonalization               : ',2f12.4,/, &
            t3,'                                   ',24('='),/,  &
            t3,'  total                           : ',2f12.4)
!
1030 format(t5,'Non-symmetric Davidson iterations (tol=',d10.2,', left =',l2,', right =',l2,'):',/, &
            t5,'------------------------------------------------------------------------------------------------',/, &
            t7,'  iter  root              eigenvalue','         rms(left)              rms(right)     max ok',/, &
            t5,'------------------------------------------------------------------------------------------------')
!
1040 format(t9,i4,2x,i4,f24.12,2d12.4,2d12.4,l3)
!
1050 format(t5,'----------------------------------------',/,&
            t7,'# target vectors:    ',i4,/,&
            t7,'# new vectors added: ',i4,/,&
            t7,'# converged vectors: ',i4,/,&
            t5,'----------------------------------------')
!


!
  end subroutine davidson_nosym_driver

  subroutine sort_eigenpairs(m,wr,wl,vr,vl,n_want,ldv,ignore,thresh,mask_in)
!
!   sort m real & imaginary eigenvalues and right & left eigenvectors of length n 
!   in decreasing order according to the real eigenvalues in the rang.e of n_want
! 
    implicit none
    integer,  intent(in)      :: m, ldv, n_want
    real(dp), intent(inout)   :: wr(m), wl(m), vr(ldv,m), vl(ldv,m)
    real(dp), intent(in)      :: thresh
    logical,  intent(in)      :: ignore
    logical,  optional        :: mask_in(m)
!   
!   local variables
!
    integer                   :: i, j, idx, min_idx(1), fin
    logical                   :: mask(m)
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
        call swap_eigenpairs(fin,idx,m,wr,wl,vr,vl,ldv)
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
      call swap_eigenpairs(i,idx,m,wr,wl,vr,vl,ldv)
! 
!
    end do
!
  end subroutine sort_eigenpairs

  subroutine swap_eigenpairs(i,j,m,wr,wl,vr,vl,ldv)
!
!   swaps m real & immaginary eigenvalues and eigenvectors of length l of the
!   indices i and j with each other 
!
    implicit none
    integer, intent(in)       :: m, ldv, i, j
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
  end subroutine swap_eigenpairs

end module mod_davidson_nosym_driver