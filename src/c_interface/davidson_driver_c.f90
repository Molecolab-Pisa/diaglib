module mod_davidson_driver_c
    use dgl_utils_c
    implicit none

contains

    subroutine davidson_driver_c(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                                 verbose, tol, max_iter, dav_iter, &
                                 shift, memory, memory_unit, metvec) bind(C, name="davidson_driver_c")
        use iso_c_binding
        implicit none

        ! C-compatible arguments
        integer(C_INT), intent(in) :: n, n_targ, n_max
        type(C_FUNPTR) :: matvec, precnd

        logical(C_BOOL), optional, intent(in) :: verbose
        integer(C_INT), optional, intent(in) :: max_iter, dav_iter
        real(C_DOUBLE), optional, intent(in) :: tol, shift
        integer(C_INT), optional, intent(in) :: memory
        type(C_PTR), optional, intent(in) :: memory_unit
        !! memory_unity has one more space to fit the line termination character in C
        type(C_FUNPTR), optional :: metvec
        
        real(C_DOUBLE), intent(inout) :: eig(n_max)
        real(C_DOUBLE), intent(inout) :: evec(n, n_max)
        logical(C_BOOL), intent(out) :: ok
        
        ! Conversione
        character(len=3) :: memory_unit_f
        logical :: verbose_f, ok_f
        ! external :: matvec_wrapper, precnd_wrapper

        ! Procedure pointer
        procedure(), pointer :: matvec_ptr, precnd_ptr

        if (present(memory_unit)) call C_string_ptr_to_F_string(memory_unit, memory_unit_f)

        ! Associa i puntatori C
        call c_f_procpointer(matvec, matvec_ptr)
        call c_f_procpointer(precnd, precnd_ptr)

        ! Conversione booleano
        verbose_f = verbose

        ! Chiamata al driver
        call dgl_davidson_driver(n, n_targ, n_max, matvec_wrapper, precnd_wrapper, eig, evec, ok_f, &
                                 dgl_verbose=verbose_f, &
                                 dgl_max_iter=max_iter, &
                                 dgl_dav_iter=dav_iter, &
                                 dgl_shift=shift, &
                                 dgl_tol=tol, &
                                 dgl_memory=memory, &
                                 dgl_memory_unit=memory_unit_f &
                                 )
        ok = ok_f

    end subroutine davidson_driver_c

    subroutine matvec_wrapper(n, m, x, ax)
        implicit none
        integer, intent(in) :: n, m
        real(dp), intent(in) :: x(n, m)
        real(dp), intent(out) :: ax(n, m)
        call matvec_ptr(n, m, x, ax)
    end subroutine

    subroutine precnd_wrapper(n, m, shift, r, z)
        implicit none
        integer, intent(in) :: n, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: r(n, m)
        real(dp), intent(out) :: z(n, m)
        call precnd_ptr(n, m, shift, r, z)
    end subroutine

  end module mod_davidson_driver_c
