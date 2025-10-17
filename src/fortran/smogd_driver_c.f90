subroutine smogd_driver_c(verbose, n, n2, n_targ, n_max, max_iter, tol, maxdav, &
                          apbmul, ambmul, spdmul, smdmul, lrprec, eig, evec, ok) bind(C, name="smogd_driver_c")
  use iso_c_binding
  use real_precision
  use diaglib
  implicit none

  ! Argomenti C
  logical(C_BOOL), value :: verbose
  integer(C_INT), value :: n, n2, n_targ, n_max, max_iter, maxdav
  real(C_DOUBLE), value :: tol
  type(C_FUNPTR), value :: apbmul, ambmul, spdmul, smdmul, lrprec
  real(C_DOUBLE), intent(inout) :: eig(n_max)
  real(C_DOUBLE), intent(inout) :: evec(n2,n_max)
  logical(C_BOOL), intent(out) :: ok

  ! Conversione
  logical :: verbose_f, ok_f
  procedure(), pointer :: apb_ptr, amb_ptr, spd_ptr, smd_ptr, prec_ptr

  ! Wrapper locali
! external :: apb_wrapper, amb_wrapper, spd_wrapper, smd_wrapper, prec_wrapper

  ! Associa i puntatori
  call c_f_procpointer(apbmul, apb_ptr)
  call c_f_procpointer(ambmul, amb_ptr)
  call c_f_procpointer(spdmul, spd_ptr)
  call c_f_procpointer(smdmul, smd_ptr)
  call c_f_procpointer(lrprec, prec_ptr)

  verbose_f = verbose

  ! Chiamata alla routine Fortran
  call smogd_driver(verbose_f, n, n2, n_targ, n_max, max_iter, tol, maxdav, &
                    apb_wrapper, amb_wrapper, spd_wrapper, smd_wrapper, prec_wrapper, eig, evec, ok_f)

  ok = ok_f

contains

  subroutine apb_wrapper(n, m, x, ax)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
    call apb_ptr(n, m, x, ax)
  end subroutine

  subroutine amb_wrapper(n, m, x, ax)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
    call amb_ptr(n, m, x, ax)
  end subroutine

  subroutine spd_wrapper(n, m, x, ax)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
    call spd_ptr(n, m, x, ax)
  end subroutine

  subroutine smd_wrapper(n, m, x, ax)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: x(n,m)
    real(dp), intent(out) :: ax(n,m)
    call smd_ptr(n, m, x, ax)
  end subroutine

  subroutine prec_wrapper(n, m, fac, xp, xm, yp, ym)
    integer, intent(in) :: n, m
    real(dp), intent(in) :: fac
    real(dp), intent(in) :: xp(n,m), xm(n,m)
    real(dp), intent(inout) :: yp(n,m), ym(n,m)
    call prec_ptr(n, m, fac, xp, xm, yp, ym)
  end subroutine

end subroutine smogd_driver_c
