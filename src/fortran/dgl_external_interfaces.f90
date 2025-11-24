module dgl_external_interfaces
    use dgl_global_utils, only: dp
    implicit none

    interface
    
        subroutine matvec_(n,m,x,y)
            import
            implicit none
            integer,                  intent(in)    :: n, m
            real(dp), dimension(n,m), intent(in)    :: x
            real(dp), dimension(n,m), intent(inout) :: y
        end subroutine matvec_

        subroutine precnd_(n,m,shift,x,y)
            import
            implicit none
            integer,                  intent(in)    :: n, m
            real(dp),                 intent(in)    :: shift
            real(dp), dimension(n,m), intent(in)    :: x
            real(dp), dimension(n,m), intent(inout) :: y
        end subroutine precnd_

        subroutine metvec_(n,m,x,y)
            import
            implicit none
            integer,                  intent(in)    :: n, m
            real(dp), dimension(n,m), intent(in)    :: x
            real(dp), dimension(n,m), intent(inout) :: y
        end subroutine metvec_

    end interface

end module