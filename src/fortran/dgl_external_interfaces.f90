module dgl_external_interfaces
    use dgl_global_utils, only: dp
    implicit none

    interface
    
        subroutine matvec_(n,m,x,y)
        !! Interface for the external routine performing the Matrix-Vector product
            import
            implicit none
            integer,                  intent(in)    :: n
            !! Lenght of the vectors to multiply
            integer,                  intent(in)    :: m
            !! Number of vectors to multiply
            real(dp), dimension(n,m), intent(in)    :: x
            !! Input vectors
            real(dp), dimension(n,m), intent(inout) :: y
            !! Output vectors
        end subroutine matvec_

        subroutine precnd_(n,m,shift,x,y)
        !! Interface for the external routine performing the preconditining of the residuals
            import
            implicit none
            integer,                  intent(in)    :: n
            !! Lenght of the vectors to multiply
            integer,                  intent(in)    :: m
            !! Number of vectors to multiply
            real(dp),                 intent(in)    :: shift
            !! Level-shifting parameter
            real(dp), dimension(n,m), intent(in)    :: x
            !! Input vectors
            real(dp), dimension(n,m), intent(inout) :: y
            !! Output vectors
        end subroutine precnd_

        subroutine metvec_(n,m,x,y)
        !! Interface for the external routine performing the Matrix-Vector product
            import
            implicit none
            integer,                  intent(in)    :: n
            !! Lenght of the vectors to multiply
            integer,                  intent(in)    :: m
            !! Number of vectors to multiply
            real(dp), dimension(n,m), intent(in)    :: x
            !! Input vectors
            real(dp), dimension(n,m), intent(inout) :: y
            !! Output vectors
        end subroutine metvec_

    end interface

end module