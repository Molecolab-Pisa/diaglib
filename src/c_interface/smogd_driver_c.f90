module mod_smogd_driver_c
    use dgl_utils_c
    implicit none
!
! The C routines of the current call are stored in module variables, as the Fortran driver
! calls them through the wrappers below. To allow calling the drivers from inside a callback,
! they are saved at the beginning of each call and restored at the end; to allow calling the
! drivers from different OpenMP threads at the same time, they are threadprivate.
!
    type(C_FUNPTR), private, save :: apb_c = C_NULL_FUNPTR, amb_c = C_NULL_FUNPTR, spd_c = C_NULL_FUNPTR, &
                                     smd_c = C_NULL_FUNPTR, prec_c = C_NULL_FUNPTR
!$omp threadprivate(apb_c, amb_c, spd_c, smd_c, prec_c)

contains

    subroutine smogd_driver_c(n2, n_targ, n_max, apbmul, ambmul, &
                              spdmul, smdmul, lrprec, eig, evec, ok, info, &
                              verbose, tol, max_iter, dav_iter, &
                              memory, memory_unit) &
        bind(C, name="dgl_smogd_driver")
        implicit none

        ! C-compatible arguments
        integer(c_ip), value, intent(in) :: n2, n_targ, n_max
        integer(c_ip), value, intent(in) :: max_iter, dav_iter
        integer(c_ip), value, intent(in) :: memory
        logical(C_BOOL), value, intent(in) :: verbose
        real(C_DOUBLE), value, intent(in) :: tol
        type(C_PTR), value, intent(in) :: memory_unit
        type(C_FUNPTR), value :: apbmul, ambmul, spdmul, smdmul, lrprec
        !
        real(C_DOUBLE), intent(inout) :: eig(n_max)
        real(C_DOUBLE), intent(inout) :: evec(n2, n_max)
        logical(C_BOOL), intent(out) :: ok
        integer(c_ip), intent(out) :: info

        character(len=2) :: memory_unit_f
!! fixed length, as required by the Fortran drivers: shorter strings
!! are blank padded, longer ones truncated, NULL gives the default unit
        logical :: verbose_f, ok_f
        integer(ip) :: info_f
        type(C_FUNPTR) :: saved(5)
!
        memory_unit_f = c_ptr_to_f_string(memory_unit)
        verbose_f = verbose
!
        ok = .false.
        info = dgl_err_input
        if (.not. pointer_ok(apbmul, "apbmul")) return
        if (.not. pointer_ok(ambmul, "ambmul")) return
        if (.not. pointer_ok(spdmul, "spdmul")) return
        if (.not. pointer_ok(smdmul, "smdmul")) return
        if (.not. pointer_ok(lrprec, "lrprec")) return
!
        saved = [apb_c, amb_c, spd_c, smd_c, prec_c]
        apb_c = apbmul
        amb_c = ambmul
        spd_c = spdmul
        smd_c = smdmul
        prec_c = lrprec
!
        call dgl_smogd_driver(n2, n_targ, n_max, apb_wrapper, amb_wrapper, spd_wrapper, smd_wrapper, prec_wrapper, &
                              eig, evec, ok_f, &
                              dgl_verbose=verbose_f, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                              dgl_tol=tol, dgl_memory=memory, dgl_memory_unit=memory_unit_f, dgl_info=info_f)
!
        apb_c = saved(1)
        amb_c = saved(2)
        spd_c = saved(3)
        smd_c = saved(4)
        prec_c = saved(5)
!
        ok = ok_f
        info = info_f

    end subroutine smogd_driver_c

    subroutine apb_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(apb_c, f)
        call f(n, m, x, ax)
    end subroutine apb_wrapper

    subroutine amb_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(amb_c, f)
        call f(n, m, x, ax)
    end subroutine amb_wrapper

    subroutine spd_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(spd_c, f)
        call f(n, m, x, ax)
    end subroutine spd_wrapper

    subroutine smd_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(smd_c, f)
        call f(n, m, x, ax)
    end subroutine smd_wrapper

    subroutine prec_wrapper(n, m, fac, xp, xm, yp, ym)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: fac
        real(dp), intent(in) :: xp(n, m), xm(n, m)
        real(dp), intent(inout) :: yp(n, m), ym(n, m)
        procedure(c_lrprec), pointer :: f
        call c_f_procpointer(prec_c, f)
        call f(n, m, fac, xp, xm, yp, ym)
    end subroutine prec_wrapper

end module mod_smogd_driver_c
