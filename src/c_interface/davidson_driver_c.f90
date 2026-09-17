module mod_davidson_driver_c
    use dgl_utils_c
    implicit none
!
! The C routines of the current call are stored in module variables, as the Fortran driver
! calls them through the wrappers below. To allow calling the drivers from inside a callback,
! they are saved at the beginning of each call and restored at the end; to allow calling the
! drivers from different OpenMP threads at the same time, they are threadprivate.
!
    type(C_FUNPTR), private, save :: matvec_c = C_NULL_FUNPTR, precnd_c = C_NULL_FUNPTR, metvec_c = C_NULL_FUNPTR
!$omp threadprivate(matvec_c, precnd_c, metvec_c)

contains

    subroutine davidson_driver_c(n, n_targ, n_max, matvec, precnd, metvec, eig, evec, ok, info, &
                                 verbose, tol, max_iter, dav_iter, shift, precnd_shift, memory, memory_unit) &
        bind(C, name="dgl_davidson_driver")
        implicit none

        ! C-compatible arguments
        integer(c_ip), value, intent(in) :: n, n_targ, n_max
        integer(c_ip), value, intent(in) :: max_iter, dav_iter
        integer(c_ip), value, intent(in) :: memory
        logical(C_BOOL), value, intent(in) :: verbose, precnd_shift
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
        logical :: verbose_f, ok_f, precnd_shift_f
        integer(ip) :: info_f
        type(C_FUNPTR) :: saved_matvec, saved_precnd, saved_metvec
        procedure(dgl_matvec), pointer :: metvec_p
!
        memory_unit_f = c_ptr_to_f_string(memory_unit)
        verbose_f = verbose
        precnd_shift_f = precnd_shift
!
        ok = .false.
        info = dgl_err_input
        if (.not. pointer_ok(matvec, "matvec")) return
        if (.not. pointer_ok(precnd, "precnd")) return
!
        saved_matvec = matvec_c
        saved_precnd = precnd_c
        saved_metvec = metvec_c
        matvec_c = matvec
        precnd_c = precnd
        metvec_c = metvec
!
        if (c_associated(metvec)) then
            metvec_p => metvec_wrapper
            call dgl_davidson_driver(n, n_targ, n_max, matvec_wrapper, precnd_wrapper, eig, evec, ok_f, &
                                     dgl_verbose=verbose_f, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                     dgl_shift=shift, dgl_tol=tol, dgl_memory=memory, &
                                     dgl_memory_unit=memory_unit_f, metvec=metvec_p, dgl_info=info_f, &
                                     dgl_precnd_shift=precnd_shift_f)
        else
            call dgl_davidson_driver(n, n_targ, n_max, matvec_wrapper, precnd_wrapper, eig, evec, ok_f, &
                                     dgl_verbose=verbose_f, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                     dgl_shift=shift, dgl_tol=tol, dgl_memory=memory, &
                                     dgl_memory_unit=memory_unit_f, dgl_info=info_f, &
                                     dgl_precnd_shift=precnd_shift_f)
        end if
!
        matvec_c = saved_matvec
        precnd_c = saved_precnd
        metvec_c = saved_metvec
!
        ok = ok_f
        info = info_f
!
    end subroutine davidson_driver_c

    subroutine matvec_wrapper(n, m, x, ax)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(matvec_c, f)
        call f(n, m, x, ax)
    end subroutine matvec_wrapper

    subroutine metvec_wrapper(n, m, x, bx)
        implicit none
        integer(ip), intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: bx(n, m)
        procedure(c_matvec), pointer :: f
        call c_f_procpointer(metvec_c, f)
        call f(n, m, x, bx)
    end subroutine metvec_wrapper

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

end module mod_davidson_driver_c
