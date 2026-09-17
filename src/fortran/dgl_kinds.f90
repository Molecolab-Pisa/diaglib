module dgl_kinds
!* Kinds of the integers and reals used in DiagLib.
    implicit none
!
#ifdef DGL_INT_KIND_8
    integer, parameter :: ip = selected_int_kind(15)
!! Global variable holding kind for integer(ip)
#elif DGL_INT_KIND_4
    integer, parameter :: ip = selected_int_kind(8)
!! Global variable holding kind for integer(ip)
#endif
    integer, parameter :: i8 = selected_int_kind(18)
!! 64-bit integer kind, always used for memory bookkeeping to avoid overflows
    integer, parameter :: dp = selected_real_kind(15)
!! Global variable holding kind for double precision
!
end module dgl_kinds
