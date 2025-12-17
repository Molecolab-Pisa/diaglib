module dgl_minor_utils
  use dgl_orthogonalizations
  implicit none
  contains

  subroutine check_guess(n,m,evec)
  !! Utility routine that checks orthogonality of the input vectors to a diagonalization driver
    implicit none
    integer,                  intent(in)    :: n
!!     
    integer,                  intent(in)    :: m
    real(dp), dimension(n,m), intent(inout) :: evec
!
    integer               :: i, j
    real(dp)              :: fac, diag_norm, out_norm, growth
    logical               :: ok
!
    real(dp), allocatable :: overlap(:,:)
!
!   check whether evec is zero.
!
    fac = dnrm2(n*m,evec,1)
    if (fac .lt. num_thresh) then
!
!     no luck. make a random guess, then orthonormalize it.
!
      call random_number(evec)
      call ortho_cd(n,m,evec,growth,ok)
    else
!
!     compute the overlap and check that the vectors are orthonormal.
!
      call mallocate(m,m,overlap)
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
      if (abs(diag_norm - one) .gt. num_thresh .or. &
          out_norm .gt. num_thresh) then
!
!       orthogonalize the guess:
!
        call ortho_cd(n,m,evec,growth,ok)
      end if
!
      call mfree(overlap)

    end if
!
    return
  end subroutine check_guess

end module dgl_minor_utils
