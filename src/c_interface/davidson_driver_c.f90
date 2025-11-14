subroutine davidson_driver_c(verbose, n, n_targ, n_max, max_iter, max_dav, tol, shift, &
                             matvec, precnd, eig, evec, ok) bind(C, name="davidson_driver_c")
  use iso_c_binding
  use diaglib_global_utils
  use diaglib
  implicit none

  ! C-compatible arguments
  logical(C_BOOL), value :: verbose
  integer(C_INT), value :: n, n_targ, n_max, max_iter, max_dav
  real(C_DOUBLE), value :: tol, shift
  type(C_FUNPTR), value :: matvec, precnd
  real(C_DOUBLE), intent(inout) :: eig(n_max)
  real(C_DOUBLE), intent(inout) :: evec(n,n_max)
  logical(C_BOOL), intent(out) :: ok

  ! Conversione
  logical :: verbose_f, ok_f
! external :: matvec_wrapper, precnd_wrapper

  ! Procedure pointer
  procedure(), pointer :: matvec_ptr, precnd_ptr

  ! Associa i puntatori C
  call c_f_procpointer(matvec, matvec_ptr)
  call c_f_procpointer(precnd, precnd_ptr)

  ! Conversione booleano
  verbose_f = verbose

  ! Chiamata al driver
  call davidson_driver(verbose_f, n, n_targ, n_max, max_iter, tol, max_dav, shift, &
                       matvec_wrapper, precnd_wrapper, eig, evec, ok_f)

  ok = ok_f

contains

  subroutine matvec_wrapper(n, m, x, ax)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
    call matvec_ptr(n, m, x, ax)
  end subroutine

  subroutine precnd_wrapper(n, m, shift, r, z)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: shift
    real(dp), intent(in) :: r(n,m)
    real(dp), intent(out) :: z(n,m)
    call precnd_ptr(n, m, shift, r, z)
  end subroutine

end subroutine davidson_driver_c
