module dgl_interface

    use dgl_global_utils, only: dgl_real => dp, dgl_int => ip, &
                                dgl_success, dgl_err_input, dgl_err_memory, &
                                dgl_err_lapack, dgl_err_ortho, dgl_err_mismatch

!
! interfaces of the user-supplied routines: matrix(metric)-vector products and preconditioners
!
    use dgl_external_interfaces, only: dgl_matvec => matvec_, dgl_precnd => precnd_, &
                                       dgl_smogd_precnd => smogd_precnd
!
    use mod_davidson_driver, only: dgl_davidson_driver => davidson_driver
    use mod_davidson_nosym_driver, only: dgl_davidson_nosym_driver => davidson_nosym_driver
    use mod_lobpcg_driver, only: dgl_lobpcg_driver => lobpcg_driver
    use mod_smogd_driver, only: dgl_smogd_driver => smogd_driver

end module
