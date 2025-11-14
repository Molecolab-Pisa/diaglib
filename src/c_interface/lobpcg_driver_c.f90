subroutine lobpcg_driver_c(verbose, gen_eig, n, n_targ, n_max, max_iter, tol, shift, &
                           matvec, precnd, bvec, eig, evec, ok) bind(C, name="lobpcg_driver_c")
  use iso_c_binding
  use diaglib_global_utils
  use diaglib
  implicit none

  ! Argomenti C
  logical(C_BOOL), value :: verbose, gen_eig
  integer(C_INT), value :: n, n_targ, n_max, max_iter
  real(C_DOUBLE), value :: tol, shift
  type(C_FUNPTR), value :: matvec, precnd, bvec
  real(C_DOUBLE), intent(inout) :: eig(n_max)
  real(C_DOUBLE), intent(inout) :: evec(n,n_max)
  logical(C_BOOL), intent(out) :: ok

  ! Conversione
  logical :: verbose_f, gen_eig_f, ok_f
  procedure(), pointer :: matvec_ptr, precnd_ptr, bvec_ptr

  ! Wrapper locali
! external :: matvec_wrapper, precnd_wrapper, bvec_wrapper

  ! Associa i puntatori
  call c_f_procpointer(matvec, matvec_ptr)
  call c_f_procpointer(precnd, precnd_ptr)
  call c_f_procpointer(bvec, bvec_ptr)

  verbose_f = verbose
  gen_eig_f = gen_eig

  ! Chiamata alla routine Fortran
  call lobpcg_driver(verbose_f, gen_eig_f, n, n_targ, n_max, max_iter, tol, shift, &
                     matvec_wrapper, precnd_wrapper, bvec_wrapper, eig, evec, ok_f)

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

  subroutine bvec_wrapper(n, m, x, bx)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: x(n,m)
    real(dp), intent(out) :: bx(n,m)
    call bvec_ptr(n, m, x, bx)
  end subroutine

end subroutine lobpcg_driver_c
