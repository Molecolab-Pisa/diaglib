module mod_lobpcg_driver_c
    use dgl_utils_c
    implicit none

    ! Procedure pointer
    procedure(), private, pointer :: matvec_ptr, precnd_ptr, metvec_ptr

contains

    subroutine lobpcg_driver_c(n, n_targ, n_max, matvec, precnd, metvec, eig, evec, ok, info, &
                               verbose, tol, max_iter, shift, memory, memory_unit) &
        bind(C, name="dgl_lobpcg_driver")
        implicit none

        ! C-compatible arguments
        integer(c_ip), value, intent(in) :: n, n_targ, n_max
        integer(c_ip), value, intent(in) :: max_iter
        integer(c_ip), value, intent(in) :: memory
        logical(C_BOOL), value, intent(in) :: verbose
        real(C_DOUBLE), value, intent(in) :: tol, shift
        type(C_PTR), value, intent(in) :: memory_unit
        type(C_FUNPTR), value :: matvec, precnd, metvec
        !
        real(C_DOUBLE), intent(inout) :: eig(n_max)
        real(C_DOUBLE), intent(inout) :: evec(n, n_max)
        logical(C_BOOL), intent(out) :: ok
        integer(c_ip), intent(out) :: info

        character(len=2) :: memory_unit_f
!! fixed length, as required by the Fortran drivers: shorter strings
!! are blank padded, longer ones truncated, NULL gives the default unit
        logical :: verbose_f, ok_f
        integer(ip) :: info_f

        memory_unit_f = c_ptr_to_f_string(memory_unit)

!       ! Associate pointers
        ok = .false.
        info = dgl_err_input
        if (.not. pointer_ok(matvec, "matvec")) return
        if (.not. pointer_ok(precnd, "precnd")) return
        call c_f_procpointer(matvec, matvec_ptr)
        call c_f_procpointer(precnd, precnd_ptr)

!       ! Bool conversion
        verbose_f = verbose

!       ! Main driver call
        if (c_associated(metvec)) then
            call c_f_procpointer(metvec, metvec_ptr)
            call dgl_lobpcg_driver(n, n_targ, n_max, matvec_wrapper, precnd_wrapper, eig, evec, ok_f, &
                                   dgl_info=info_f, &
                                   dgl_verbose=verbose_f, &
                                   dgl_max_iter=max_iter, &
                                   dgl_shift=shift, &
                                   dgl_tol=tol, &
                                   dgl_memory=memory, &
                                   dgl_memory_unit=memory_unit_f, &
                                   metvec=metvec_ptr &
                                   )
        else
            call dgl_lobpcg_driver(n, n_targ, n_max, matvec_wrapper, precnd_wrapper, eig, evec, ok_f, &
                                   dgl_info=info_f, &
                                   dgl_verbose=verbose_f, &
                                   dgl_max_iter=max_iter, &
                                   dgl_shift=shift, &
                                   dgl_tol=tol, &
                                   dgl_memory=memory, &
                                   dgl_memory_unit=memory_unit_f &
                                   )
        end if
!       ! Chiamata al driver
        ok = ok_f
        info = info_f

    end subroutine lobpcg_driver_c

    subroutine matvec_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call matvec_ptr(n, m, x, ax)
    end subroutine

    subroutine precnd_wrapper(n, m, shift, r, z)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: r(n, m)
        real(dp), intent(inout) :: z(n, m)
        call precnd_ptr(n, m, shift, r, z)
    end subroutine

end module mod_lobpcg_driver_c
