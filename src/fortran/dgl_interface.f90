module dgl_interface

    use dgl_global_utils, only: dgl_real => dp, dgl_int => ip

    use mod_davidson_driver, only: dgl_davidson_driver => davidson_driver
    use mod_davidson_nosym_driver, only: dgl_davidson_nosym_driver => davidson_nosym_driver
    use mod_lobpcg_driver, only: dgl_lobpcg_driver => lobpcg_driver
    use mod_smogd_driver, only: dgl_smogd_driver => smogd_driver
    use dgl_orthogonalizations, only: dgl_ortho => ortho, dgl_b_ortho => b_ortho, &
                                      dgl_ortho_cd => ortho_cd, dgl_ortho_vs_x => ortho_vs_x, &
                                      dgl_b_ortho_vs_x => b_ortho_vs_x

end module
