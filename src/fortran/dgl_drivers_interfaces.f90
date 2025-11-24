module dgl_drivers_interfaces
    use dgl_global_utils, only: dp
    use dgl_external_interfaces
    implicit none

    interface

        subroutine davidson_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav,&
                                shift,matvec,precnd,eig,evec,ok,metvec)
            import
            implicit none

            logical,                      intent(in)    :: verbose
            integer,                      intent(in)    :: n, n_targ, n_max
            integer,                      intent(in)    :: max_iter, max_dav
            real(dp),                     intent(in)    :: tol, shift
            real(dp), dimension(n_max),   intent(inout) :: eig
            real(dp), dimension(n,n_max), intent(inout) :: evec
            logical,                      intent(inout) :: ok
            procedure(matvec_) :: matvec
            procedure(precnd_) :: precnd
            procedure(metvec_), pointer, optional :: metvec

        end subroutine

        subroutine lobpcg_driver(verbose,n,n_targ,n_max,max_iter,tol, &
                           shift,matvec,precnd,eig,evec,ok,metvec)
        import
        implicit none
            logical,                      intent(in)    :: verbose
            integer,                      intent(in)    :: n, n_targ, n_max
            integer,                      intent(in)    :: max_iter
            real(dp),                     intent(in)    :: tol, shift
            real(dp), dimension(n_max),   intent(inout) :: eig
            real(dp), dimension(n,n_max), intent(inout) :: evec
            logical,                      intent(inout) :: ok
            procedure(matvec_) :: matvec
            procedure(precnd_) :: precnd
            procedure(metvec_), pointer, optional :: metvec

        end subroutine lobpcg_driver

    end interface

end module dgl_drivers_interfaces