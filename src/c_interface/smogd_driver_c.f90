module mod_smogd_driver_c
    use dgl_utils_c
    implicit none

    ! Procedure pointer
    procedure(), private, pointer :: apb_ptr, amb_ptr, spd_ptr, smd_ptr, prec_ptr

contains

    subroutine smogd_driver_c(n2, n_targ, n_max, apbmul, ambmul, &
                              spdmul, smdmul, lrprec, eig, evec, ok, &
                              verbose, tol, max_iter, dav_iter, &
                              memory, memory_unit) &
        bind(C, name="dgl_smogd_driver")
        implicit none

        ! C-compatible arguments
#ifdef DGL_INT_KIND_4
        integer(C_INT), value, intent(in) :: n2, n_targ, n_max
        integer(C_INT), value, intent(in) :: max_iter, dav_iter
        integer(C_INT), value, intent(in) :: memory
#elif DGL_INT_KIND_8
        integer(C_LONG), value, intent(in) :: n2, n_targ, n_max
        integer(C_LONG), value, intent(in) :: max_iter, dav_iter
        integer(C_LONG), value, intent(in) :: memory
#endif
        logical(C_BOOL), value, intent(in) :: verbose
        real(C_DOUBLE), value, intent(in) :: tol
        type(C_PTR), value, intent(in) :: memory_unit
        type(C_FUNPTR), value :: apbmul, ambmul, spdmul, smdmul, lrprec
        !
        real(C_DOUBLE), intent(inout) :: eig(n_max)
        real(C_DOUBLE), intent(inout) :: evec(n2, n_max)
        logical(C_BOOL), intent(out) :: ok

        character(len=:), allocatable :: memory_unit_f
        logical :: verbose_f, ok_f

        memory_unit_f = c_ptr_to_f_string(memory_unit)

        ! Associa i puntatori
        call check_pointer(apbmul, "apbmul")
        call check_pointer(ambmul, "ambmul")
        call check_pointer(spdmul, "spdmul")
        call check_pointer(smdmul, "smdmul")
        call check_pointer(lrprec, "lrprec")
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
                              dgl_verbose=verbose_f, &
                              dgl_max_iter=max_iter, &
                              dgl_dav_iter=dav_iter, &
                              dgl_tol=tol, &
                              dgl_memory=memory, &
                              dgl_memory_unit=memory_unit_f &
                              )

        ok = ok_f

    end subroutine smogd_driver_c

    subroutine apb_wrapper(n, m, x, ax)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call apb_ptr(n, m, x, ax)
    end subroutine

    subroutine amb_wrapper(n, m, x, ax)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call amb_ptr(n, m, x, ax)
    end subroutine

    subroutine spd_wrapper(n, m, x, ax)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call spd_ptr(n, m, x, ax)
    end subroutine

    subroutine smd_wrapper(n, m, x, ax)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call smd_ptr(n, m, x, ax)
    end subroutine

    subroutine prec_wrapper(n, m, fac, xp, xm, yp, ym)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: fac
        real(dp), intent(in) :: xp(n, m), xm(n, m)
        real(dp), intent(inout) :: yp(n, m), ym(n, m)
        call prec_ptr(n, m, fac, xp, xm, yp, ym)
    end subroutine

end module mod_smogd_driver_c
