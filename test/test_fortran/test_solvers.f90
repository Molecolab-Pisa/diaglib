module solvers
use dgl_interface
use matvecs

  integer :: i, j
  logical :: ok
  integer,  parameter :: lutest=100
  real(dgl_real), allocatable :: eig(:), evec(:,:), evec_l(:,:)
  integer,  parameter :: memory = 100
  character(len=2),parameter :: memory_unit = "MB"
  procedure(), pointer :: mx_p => null()


contains

subroutine test_davidson(n, n_targ, n_max, verbose, max_iter, dav_iter, tol, shift)
  implicit none
  integer, intent(in) :: n, n_targ, n_max
  integer, intent(in), optional :: max_iter, dav_iter
  real(dgl_real), intent(in), optional :: tol, shift
  logical, intent(in), optional :: verbose
!
!
! open a text file for the output, to be used to compare the results with a reference.
!
  open (file='output_fortran.txt',form='formatted',access='sequential', & 
        unit=lutest,status='unknown')
!
!
  allocate (eig(n_max), evec(n, n_max))
  ok = .false.
!
! test davidson:
!
  eig  = dgl_zero
  evec = dgl_zero
  do i = 1, n_max
    evec(i,i) = dgl_one
  end do
  ok = .false.
!
  write(6,*) ' testing Davidson:'
  call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                        dgl_verbose = verbose, &
                        dgl_max_iter = max_iter, &
                        dgl_dav_iter = dav_iter, &
                        dgl_tol = tol, &
                        dgl_shift = shift, &
                        dgl_memory = memory, &
                        dgl_memory_unit = memory_unit)
!
  if (ok) then
    write(6,*) ' Davidson converged.'
    write(lutest,1000) 'Davidson'
    write(lutest,*)
    write(lutest,1010)
    do i = 1, n_targ
      write(lutest,1020) i, eig(i)
    end do
    do i = 1, n_targ
      if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
    end do
    write(lutest,1031) (i, i = 1, n_targ)
    do j = 1, n
      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
    end do
    write(lutest,*)
  else
    write(6,*) ' Davidson failed to converge.'
  end if
!
  deallocate(eig,evec)

 
end subroutine test_davidson


end module solvers