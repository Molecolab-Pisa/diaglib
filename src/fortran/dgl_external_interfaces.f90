module dgl_external_interfaces
!* Module contaning all the interfaces for the external procedure,
! like matrix-vector multipltications.
    use dgl_global_utils, only: dp
    implicit none

    interface
    !! Interface to external procedures used in diaglib
        subroutine matvec_(n, m, x, y)
        !! Interface for the external routine performing the Matrix-Vector product
            import
            implicit none
            integer, intent(in) :: n
            !! Lenght of the vectors to multiply
            integer, intent(in) :: m
            !! Number of vectors to multiply
            real(dp), dimension(n, m), intent(in) :: x
            !! Input vectors
            real(dp), dimension(n, m), intent(inout) :: y
            !! Output vectors
        end subroutine matvec_

        subroutine precnd_(n, m, shift, x, y)
        !! Interface for the external routine performing the preconditining of the residuals
            import
            implicit none
            integer, intent(in) :: n
            !! Lenght of the vectors to multiply
            integer, intent(in) :: m
            !! Number of vectors to multiply
            real(dp), intent(in) :: shift
            !! Level-shifting parameter
            real(dp), dimension(n, m), intent(in) :: x
            !! Input vectors
            real(dp), dimension(n, m), intent(inout) :: y
            !! Output vectors
        end subroutine precnd_

        subroutine metvec_(n, m, x, y)
        !! Interface for the external routine performing the Matrix-Vector product
            import
            implicit none
            integer, intent(in) :: n
            !! Lenght of the vectors to multiply
            integer, intent(in) :: m
            !! Number of vectors to multiply
            real(dp), dimension(n, m), intent(in) :: x
            !! Input vectors
            real(dp), dimension(n, m), intent(inout) :: y
            !! Output vectors
        end subroutine metvec_

        subroutine smogd_matvec(n, m, x, y)
        !! Interface for the external routine performing the Matrix-Vector products for SMO-GD
            import
            implicit none
            integer, intent(in) :: n
            !! Lenght of the vectors to multiply
            integer, intent(in) :: m
            !! Number of vectors to multiply
            real(dp), dimension(n, m), intent(in) :: x
            !! Input vectors
            real(dp), dimension(n, m), intent(inout) :: y
            !! Output vectors
        end subroutine smogd_matvec

        subroutine smogd_precnd(n, m, fac, xp, xm, yp, ym)
            import
            implicit none
            integer, intent(in) :: n
            !! Lenght of the vectors to multiply
            integer, intent(in) :: m
            !! Number of vectors to multiply
            real(dp), intent(in) :: fac
            !! ???
            real(dp), dimension(n, m), intent(in) :: xp
            !! Input vectors plus vectors
            real(dp), dimension(n, m), intent(in) :: xm
            !! Input vectors minus vectors
            real(dp), dimension(n, m), intent(inout) :: yp
            !! Output vectors plus vectors
            real(dp), dimension(n, m), intent(inout) :: ym
            !! Output vectors minus vectors
        end subroutine smogd_precnd

    end interface

end module
