module mod_davidson_nosym_driver_c
    use dgl_utils_c
    implicit none

    ! Procedure pointer
    procedure(), private, pointer :: matvec_r_ptr, matvec_l_ptr, precnd_ptr

contains

    subroutine davidson_nosym_driver_c(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, &
                                       eig, evec_1, evec_2, ok, info, &
                                       verbose, tol, max_iter, dav_iter, &
                                       shift, memory, memory_unit) &
        bind(C, name="dgl_davidson_nosym_driver")
        implicit none

        ! C-compatible arguments
        integer(c_ip), value, intent(in) :: n, n_targ, n_max
        integer(c_ip), value, intent(in) :: max_iter, dav_iter
        integer(c_ip), value, intent(in) :: memory
        logical(C_BOOL), value, intent(in) :: verbose
        real(C_DOUBLE), value, intent(in) :: tol, shift
        type(C_PTR), value, intent(in) :: memory_unit, side
        type(C_PTR), value, intent(in) :: evec_2
!! only used (and required) if side = "LR", may be NULL otherwise
        type(C_FUNPTR), value :: matvec_r, matvec_l, precnd
        !
        real(C_DOUBLE), intent(inout) :: eig(n_max)
        real(C_DOUBLE), intent(inout) :: evec_1(n, n_max)
        logical(C_BOOL), intent(out) :: ok
        integer(c_ip), intent(out) :: info

        character(len=2) :: memory_unit_f, side_f
!! fixed length, as required by the Fortran driver: shorter strings
!! are blank padded, longer ones truncated, NULL gives the default unit
        logical :: verbose_f, ok_f
        integer(ip) :: info_f
        real(C_DOUBLE), pointer :: evec_2_f(:, :)
!
        memory_unit_f = c_ptr_to_f_string(memory_unit)
        side_f = c_ptr_to_f_string(side)
!
!       ! Associate pointers
        ok = .false.
        info = dgl_err_input
        if (.not. pointer_ok(matvec_r, "matvec_r")) return
        if (.not. pointer_ok(matvec_l, "matvec_l")) return
        if (.not. pointer_ok(precnd, "precnd")) return
        call c_f_procpointer(matvec_r, matvec_r_ptr)
        call c_f_procpointer(matvec_l, matvec_l_ptr)
        call c_f_procpointer(precnd, precnd_ptr)

        verbose_f = verbose

        ! Chiamata al driver
        if (trim(side_f) == "LR") then
            if (.not. pointer_ok(evec_2, "evec_2")) return
            call c_f_pointer(evec_2, evec_2_f, [n, n_max])
            call dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r_wrapper, matvec_l_wrapper, precnd_wrapper, &
                                           side_f, eig, evec_1, ok_f, &
                                           evec_2=evec_2_f, &
                                           dgl_info=info_f, &
                                           dgl_verbose=verbose_f, &
                                           dgl_max_iter=max_iter, &
                                           dgl_dav_iter=dav_iter, &
                                           dgl_shift=shift, &
                                           dgl_tol=tol, &
                                           dgl_memory=memory, &
                                           dgl_memory_unit=memory_unit_f &
                                           )
        else
            call dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r_wrapper, matvec_l_wrapper, precnd_wrapper, &
                                           side_f, eig, evec_1, ok_f, &
                                           dgl_info=info_f, &
                                           dgl_verbose=verbose_f, &
                                           dgl_max_iter=max_iter, &
                                           dgl_dav_iter=dav_iter, &
                                           dgl_shift=shift, &
                                           dgl_tol=tol, &
                                           dgl_memory=memory, &
                                           dgl_memory_unit=memory_unit_f &
                                           )
        end if

        ok = ok_f
        info = info_f

    end subroutine davidson_nosym_driver_c

    subroutine matvec_l_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
!
        call matvec_l_ptr(n, m, x, ax)
!
    end subroutine

    subroutine matvec_r_wrapper(n, m, x, ax)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
!
        call matvec_r_ptr(n, m, x, ax)
!
    end subroutine

    subroutine precnd_wrapper(n, m, shift, r, z)
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: r(n, m)
        real(dp), intent(inout) :: z(n, m)
!
        call precnd_ptr(n, m, shift, r, z)
!
    end subroutine

end module mod_davidson_nosym_driver_c
