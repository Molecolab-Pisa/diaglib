module dgl_drivers_interfaces
    use dgl_global_utils, only: dp
    use dgl_external_interfaces
    implicit none

    interface

        subroutine davidson_driver(n,n_targ,n_max,matvec,precnd,eig,evec,ok, &
              dgl_verbose, dgl_max_iter, dgl_tol, dgl_max_dav, &
              dgl_shift, dgl_memory, metvec)
            import
            implicit none
            integer,                      intent(in)    :: n, n_targ, n_max
            real(dp), dimension(n_max),   intent(inout) :: eig
            real(dp), dimension(n,n_max), intent(inout) :: evec
            logical,                      intent(inout) :: ok
            procedure(matvec_) :: matvec
            procedure(precnd_) :: precnd

            logical,  optional,            intent(in)    :: dgl_verbose
            integer,  optional,            intent(in)    :: dgl_max_iter, dgl_max_dav, dgl_memory
            real(dp), optional,            intent(in)    :: dgl_tol, dgl_shift
            procedure(metvec_), pointer, optional :: metvec

        end subroutine davidson_driver

        subroutine lobpcg_driver(n,n_targ,n_max,matvec,precnd,eig,evec,ok, &
              dgl_verbose, dgl_max_iter, dgl_tol, &
              dgl_shift, dgl_memory, metvec)
            import
            implicit none
            integer,                      intent(in)    :: n, n_targ, n_max
            real(dp), dimension(n_max),   intent(inout) :: eig
            real(dp), dimension(n,n_max), intent(inout) :: evec
            logical,                      intent(inout) :: ok
            procedure(matvec_) :: matvec
            procedure(precnd_) :: precnd

            logical,  optional,            intent(in)    :: dgl_verbose
            integer,  optional,            intent(in)    :: dgl_max_iter, dgl_memory
            real(dp), optional,            intent(in)    :: dgl_tol, dgl_shift
            procedure(metvec_), pointer, optional :: metvec

        end subroutine lobpcg_driver
        
    end interface

end module dgl_drivers_interfaces