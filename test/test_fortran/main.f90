program main
    use dgl_interface
    use direct_solvers
    use direct_matvecs
    use sparsematrices

    implicit none

    !call reset_output_file()

    !call test_davidson(500, 10, 10, verbose=.false.)
    !call test_lobpcg(500, 10, 10, verbose=.false.)
    !call test_nosym_davidson(500, 10, 10, "R", verbose=.true.)
    !call test_smogd(500, 10, 10, verbose=.false.)

    call read_dimensions(trim(matrix_file(1)))
    call read_matrix(trim(matrix_file(1)))

end program main
