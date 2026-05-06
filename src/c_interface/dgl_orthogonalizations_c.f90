module dgl_orthogonalizations_c
    use dgl_utils_c
    implicit none

contains

    subroutine ortho_cd_c(n, m, u, growth, ok) bind(C,name="dgl_ortho_cd")
        implicit none
!        
#ifdef DGL_INT_KIND_4
        integer(C_INT), value, intent(in) :: n
        integer(C_INT), value, intent(in) :: m
#elif DGL_INT_KIND_8
        integer(C_LONG), value, intent(in) :: n
        integer(C_LONG), value, intent(in) :: m
#endif DGL_INT_KIND_4
        real(C_DOUBLE), dimension(n, m), intent(inout) :: u
        real(C_DOUBLE), intent(out) :: growth
        logical(C_BOOL), intent(out) :: ok

        logical :: ok_f
        
        call dgl_ortho_cd(n, m, u, growth, ok_f)

        ok = ok_f

    end subroutine ortho_cd_c

    subroutine ortho_vs_x_c(n, m, k, x, u) bind(C,name="dgl_ortho_vs_x")
        implicit none

#ifdef DGL_INT_KIND_4
        integer(C_INT), value, intent(in) :: n
        integer(C_INT), value, intent(in) :: m
        integer(C_INT), value, intent(in) :: k
#elif DGL_INT_KIND_8
        integer(C_LONG), value, intent(in) :: n
        integer(C_LONG), value, intent(in) :: m
        integer(C_LONG), value, intent(in) :: k
#endif DGL_INT_KIND_4

        real(C_DOUBLE), dimension(n, m), intent(in) :: x
        real(C_DOUBLE), dimension(n, k), intent(inout) :: u

        call dgl_ortho_vs_x(n, m, k, x, u)

    end subroutine ortho_vs_x_c

    subroutine b_ortho_vs_x_c(n, m, k, x, bx, u) bind(C,name="dgl_b_ortho_vs_x")
        implicit none
        
#ifdef DGL_INT_KIND_4
        integer(C_INT), value, intent(in) :: n
        integer(C_INT), value, intent(in) :: m
        integer(C_INT), value, intent(in) :: k
#elif DGL_INT_KIND_8
        integer(C_LONG), value, intent(in) :: n
        integer(C_LONG), value, intent(in) :: m
        integer(C_LONG), value, intent(in) :: k
#endif DGL_INT_KIND_4

        real(C_DOUBLE), dimension(n, m), intent(in) :: x
        real(C_DOUBLE), dimension(n, m), intent(in) :: bx
        real(C_DOUBLE), dimension(n, k), intent(inout) :: u

        call dgl_b_ortho_vs_x(n, m, k, x, bx, u)

    end subroutine b_ortho_vs_x_c

end module dgl_orthogonalizations_c