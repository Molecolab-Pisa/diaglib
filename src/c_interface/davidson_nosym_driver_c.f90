module mod_davidson_nosym_driver_c
    use dgl_utils_c
    implicit none
!
! The C routines of the current call are stored in module variables, as the Fortran driver
! calls them through the wrappers below. To allow calling the drivers from inside a callback,
! they are saved at the beginning of each call and restored at the end; to allow calling the
! drivers from different OpenMP threads at the same time, they are threadprivate.
!
    type(C_FUNPTR), private, save :: matvec_r_c = C_NULL_FUNPTR, matvec_l_c = C_NULL_FUNPTR, precnd_c = C_NULL_FUNPTR
!$omp threadprivate(matvec_r_c, matvec_l_c, precnd_c)

contains

    subroutine davidson_nosym_driver_c(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, &
                                       eig, evec_1, evec_2, ok, info, &
                                       verbose, tol, max_iter, dav_iter, &
                                       shift, precnd_shift, memory, memory_unit) &
        bind(C, name="dgl_davidson_nosym_driver")
        implicit none

        ! C-compatible arguments
        integer(c_ip), value, intent(in) :: n, n_targ, n_max
        integer(c_ip), value, intent(in) :: max_iter, dav_iter
        integer(c_ip), value, intent(in) :: memory
        logical(C_BOOL), value, intent(in) :: verbose, precnd_shift
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
        logical :: verbose_f, ok_f, precnd_shift_f
        integer(ip) :: info_f
        real(C_DOUBLE), pointer :: evec_2_f(:, :)
        type(C_FUNPTR) :: saved_matvec_r, saved_matvec_l, saved_precnd
!
        memory_unit_f = c_ptr_to_f_string(memory_unit)
        side_f = c_ptr_to_f_string(side)
        verbose_f = verbose
        precnd_shift_f = precnd_shift
!
        ok = .false.
        info = dgl_err_input
        if (.not. pointer_ok(matvec_r, "matvec_r")) return
        if (.not. pointer_ok(matvec_l, "matvec_l")) return
        if (.not. pointer_ok(precnd, "precnd")) return
        if (trim(side_f) == "LR") then
            if (.not. pointer_ok(evec_2, "evec_2")) return
        end if
!
        saved_matvec_r = matvec_r_c
        saved_matvec_l = matvec_l_c
        saved_precnd = precnd_c
        matvec_r_c = matvec_r
        matvec_l_c = matvec_l
        precnd_c = precnd
!
        if (trim(side_f) == "LR") then
            call c_f_pointer(evec_2, evec_2_f, [n, n_max])
            call dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r_wrapper, matvec_l_wrapper, precnd_wrapper, &
                                           side_f, eig, evec_1, ok_f, evec_2=evec_2_f, &
                                           dgl_verbose=verbose_f, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                           dgl_shift=shift, dgl_tol=tol, dgl_memory=memory, &
                                           dgl_memory_unit=memory_unit_f, dgl_info=info_f, &
                                           dgl_precnd_shift=precnd_shift_f)
        else
            call dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r_wrapper, matvec_l_wrapper, precnd_wrapper, &
                                           side_f, eig, evec_1, ok_f, &
                                           dgl_verbose=verbose_f, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                           dgl_shift=shift, dgl_tol=tol, dgl_memory=memory, &
                                           dgl_memory_unit=memory_unit_f, dgl_info=info_f, &
                                           dgl_precnd_shift=precnd_shift_f)
        end if
!
        matvec_r_c = saved_matvec_r
        matvec_l_c = saved_matvec_l
        precnd_c = saved_precnd
!
        ok = ok_f
        info = info_f

    end subroutine davidson_nosym_driver_c

    subroutine matvec_r_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(matvec_r_c, f)
        call f(n, m, x, ax)
    end subroutine matvec_r_wrapper

    subroutine matvec_l_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(matvec_l_c, f)
        call f(n, m, x, ax)
    end subroutine matvec_l_wrapper

    subroutine precnd_wrapper(n, m, shift, r, z)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: r(n, m)
        real(dp), intent(inout) :: z(n, m)
        procedure(c_precnd), pointer :: f
        call c_f_procpointer(precnd_c, f)
        call f(n, m, shift, r, z)
    end subroutine precnd_wrapper

end module mod_davidson_nosym_driver_c
