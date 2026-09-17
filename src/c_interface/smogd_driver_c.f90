module mod_smogd_driver_c
    use dgl_utils_c
    implicit none

    ! Procedure pointer
    procedure(), private, pointer :: apb_ptr, amb_ptr, spd_ptr, smd_ptr, prec_ptr

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

        memory_unit_f = c_ptr_to_f_string(memory_unit)

        ! Associa i puntatori
        ok = .false.
        info = dgl_err_input
        if (.not. pointer_ok(apbmul, "apbmul")) return
        if (.not. pointer_ok(ambmul, "ambmul")) return
        if (.not. pointer_ok(spdmul, "spdmul")) return
        if (.not. pointer_ok(smdmul, "smdmul")) return
        if (.not. pointer_ok(lrprec, "lrprec")) return
        call c_f_procpointer(apbmul, apb_ptr)
        call c_f_procpointer(ambmul, amb_ptr)
        call c_f_procpointer(spdmul, spd_ptr)
        call c_f_procpointer(smdmul, smd_ptr)
        call c_f_procpointer(lrprec, prec_ptr)

        verbose_f = verbose

        ! Call to the driver
        call dgl_smogd_driver(n2, n_targ, n_max, apb_wrapper, amb_wrapper, &
                              spd_wrapper, smd_wrapper, prec_wrapper, &
                              eig, evec, ok_f, &
                              dgl_info=info_f, &
                              dgl_verbose=verbose_f, &
                              dgl_max_iter=max_iter, &
                              dgl_dav_iter=dav_iter, &
                              dgl_tol=tol, &
                              dgl_memory=memory, &
                              dgl_memory_unit=memory_unit_f &
                              )

        ok = ok_f
        info = info_f

    end subroutine smogd_driver_c

    subroutine apb_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call apb_ptr(n, m, x, ax)
    end subroutine

    subroutine amb_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call amb_ptr(n, m, x, ax)
    end subroutine

    subroutine spd_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call spd_ptr(n, m, x, ax)
    end subroutine

    subroutine smd_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call smd_ptr(n, m, x, ax)
    end subroutine

    subroutine prec_wrapper(n, m, fac, xp, xm, yp, ym)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: fac
        real(dp), intent(in) :: xp(n, m), xm(n, m)
        real(dp), intent(inout) :: yp(n, m), ym(n, m)
        call prec_ptr(n, m, fac, xp, xm, yp, ym)
    end subroutine

end module mod_smogd_driver_c
