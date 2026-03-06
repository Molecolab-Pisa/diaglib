program main
    use dgl_interface
    use direct_solvers
    use direct_matvecs

    implicit none

    call reset_output_file()

    call test_davidson(500, 5, 10, verbose=.false., tol = 1.d-12)
    call test_lobpcg(500, 5, 10, verbose=.false., tol = 1.d-12)
    call test_nosym_davidson(500, 5, 10, "LR", verbose=.true., tol = 1.d-12)
    call test_smogd(500, 5, 10, verbose=.false., tol = 1.d-12)


end program main