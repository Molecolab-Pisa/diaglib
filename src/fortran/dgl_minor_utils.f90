module dgl_minor_utils
  use dgl_orthogonalizations
  implicit none
  contains

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

end module dgl_minor_utils
