module mod_davidson_driver_c
    use dgl_utils_c
    implicit none

    ! Procedure pointer
    procedure(), private, pointer :: matvec_ptr, precnd_ptr, metvec_ptr

contains

    subroutine davidson_driver_c(n, n_targ, n_max, matvec, precnd, metvec, eig, evec, ok, &
                                 verbose, tol, max_iter, dav_iter, shift, memory, memory_unit) &
                                 bind(C, name="davidson_driver_c")
        implicit none

        ! C-compatible arguments
#ifdef DGL_INT_KIND_4
        integer(C_INT), value, intent(in) :: n, n_targ, n_max
        integer(C_INT), value, intent(in) :: max_iter, dav_iter
        integer(C_INT), value, intent(in) :: memory
#elif DGL_INT_KIND_8
        integer(C_LONG), value, intent(in) :: n, n_targ, n_max
        integer(C_LONG), value, intent(in) :: max_iter, dav_iter
        integer(C_LONG), value, intent(in) :: memory
#endif
        logical(C_BOOL), value, intent(in) :: verbose
        real(C_DOUBLE), value, intent(in) :: tol, shift
        type(C_PTR), value, intent(in) :: memory_unit
        type(C_FUNPTR), value :: matvec, precnd, metvec
        !
        real(C_DOUBLE), intent(inout) :: eig(n_max)
        real(C_DOUBLE), intent(inout) :: evec(n, n_max)
        logical(C_BOOL), intent(out) :: ok

        character(len=:), allocatable :: memory_unit_f
        logical :: verbose_f, ok_f
!
        memory_unit_f = c_ptr_to_f_string(memory_unit)
!
!       ! Associate pointers
        call c_f_procpointer(matvec, matvec_ptr)
        call c_f_procpointer(precnd, precnd_ptr)
        call c_f_procpointer(metvec, metvec_ptr)
!
!       ! Bool conversion
        verbose_f = verbose
!
!       ! Main driver call
        call dgl_davidson_driver(n, n_targ, n_max, matvec_wrapper, precnd_wrapper, eig, evec, ok_f, &
                                 dgl_verbose=verbose_f, &
                                 dgl_max_iter=max_iter, &
                                 dgl_dav_iter=dav_iter, &
                                 dgl_shift=shift, &
                                 dgl_tol=tol, &
                                 dgl_memory=memory, &
                                 dgl_memory_unit=memory_unit_f, &
                                 metvec=metvec_ptr &
                                 )
        ! Bool conversion
        ok = ok_f
!
    end subroutine davidson_driver_c

    subroutine matvec_wrapper(n, m, x, ax)
        implicit none
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(inout) :: ax(n, m)
        call matvec_ptr(n, m, x, ax)
    end subroutine

    subroutine precnd_wrapper(n, m, shift, r, z)
        implicit none
        integer, intent(in) :: n, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: r(n, m)
        real(dp), intent(inout) :: z(n, m)
        call precnd_ptr(n, m, shift, r, z)
    end subroutine

end module mod_davidson_driver_c
