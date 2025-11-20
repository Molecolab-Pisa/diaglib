module dgl_interfaces
  implicit none
  interface
        subroutine myproc(a)
            real :: a
        end subroutine myproc

        subroutine davidson_driver(verbose,n,n_targ,n_max,max_iter,tol,max_dav,&
                                shift,matvec,precnd,eig,evec,ok,proc_pointer)
            use dgl_minor_utils
            implicit none

            logical,                      intent(in)    :: verbose
            integer,                      intent(in)    :: n, n_targ, n_max
            integer,                      intent(in)    :: max_iter, max_dav
            real(dp),                     intent(in)    :: tol, shift
            real(dp), dimension(n_max),   intent(inout) :: eig
            real(dp), dimension(n,n_max), intent(inout) :: evec
            logical,                      intent(inout) :: ok
            external                                    :: matvec, precnd
            procedure(myproc), pointer, optional :: proc_pointer


        end subroutine
    end interface

end module dgl_interfaces