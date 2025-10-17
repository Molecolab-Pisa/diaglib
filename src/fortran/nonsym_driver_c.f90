subroutine nonsym_driver_c(verbose, n, n_targ, n_max, max_iter, tol, max_dav, shift,  &
                           matvec_r, matvec_l, precnd, eig, evec_r, evec_l, side, ok) &
                           bind(C, name="nonsym_driver_c")
  use iso_c_binding
  use real_precision
  use diaglib
  implicit none

  ! C-compatible arguments
  logical(C_BOOL), value         :: verbose
  integer(C_INT),  value         :: n, n_targ, n_max, max_iter, max_dav, side
  real(C_DOUBLE),  value         :: tol, shift
  type(C_FUNPTR),  value         :: matvec_r, matvec_l, precnd
  real(C_DOUBLE),  intent(inout) :: eig(n_max)
  real(C_DOUBLE),  intent(inout) :: evec_r(n,n_max), evec_l(n,n_max)
  logical(C_BOOL), intent(out)   :: ok

  ! Conversione
  logical :: verbose_f, ok_f

  ! Procedure pointer
  procedure(), pointer :: matvec_l_ptr, matvec_r_ptr, precnd_ptr

  ! Associa i puntatori C
  call c_f_procpointer(matvec_r, matvec_r_ptr)
  call c_f_procpointer(matvec_l, matvec_l_ptr)
  call c_f_procpointer(precnd, precnd_ptr)

  ! Conversione booleano
  verbose_f = verbose

  ! Chiamata al driver
  call nonsym_driver(verbose_f, n, n_targ, n_max, max_iter, tol, max_dav, shift, &
                     matvec_r_wrapper, matvec_l_wrapper, precnd_wrapper, eig,  &
                     evec_r, evec_l, side, ok_f)

  ok = ok_f

contains

  subroutine matvec_l_wrapper(n, m, x, ax)
    integer,  intent(in)  :: n, m
    real(dp), intent(in)  :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
!
    call matvec_l_ptr(n, m, x, ax)
!
  end subroutine

  subroutine matvec_r_wrapper(n, m, x, ax)
    integer,  intent(in)  :: n, m
    real(dp), intent(in)  :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
!
    call matvec_r_ptr(n, m, x, ax)
!
  end subroutine

  subroutine precnd_wrapper(n, m, shift, r, z)
    integer,  intent(in)  :: n, m
    real(dp), intent(in)  :: shift
    real(dp), intent(in)  :: r(n,m)
    real(dp), intent(out) :: z(n,m)
!
    call precnd_ptr(n, m, shift, r, z)
!
  end subroutine

end subroutine nonsym_driver_c
