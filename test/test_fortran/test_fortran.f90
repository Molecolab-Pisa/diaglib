program test_fortran
!!
!! tests the various functionalities of diaglib.
!!
    use dgl_interface
    use direct_matvecs
    use utility
    implicit none
    integer, parameter :: n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 10
    logical :: verbose = .false.
    real(dp), parameter :: tol = 1.0e-10_dp, shift = 1.0_dp
    integer, parameter :: memory = 100
    character(len=2), parameter :: memory_unit = "MB"
    procedure(), pointer :: mx_p => null()
!
    logical :: ok
    real(dp), allocatable :: eig(:), evec(:, :), evec_2(:, :)
    integer, parameter :: lutest = 100
    character(len=30) :: results_fname = "output_fortran.txt"
!
! open a text file for the output, to be used to compare the results with a reference.
!
    open (file=trim(results_fname), form='formatted', access='sequential', &
          unit=lutest, status='unknown')
!
! allocate space for eigenvalues and eigenvectors
!
    allocate (eig(n_max), evec(n, n_max), evec_2(n, n_max))
!
! Associate pointer for generalized problems to metric-vector
! product routine
!
    mx_p => mx
!
! test davidson:
!
    call init_eigenpairs(n, n_max, eig, evec, ok)
!
    write (6, *) ' testing Davidson:'
    call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                             dgl_verbose=verbose, &
                             dgl_max_iter=max_iter, &
                             dgl_dav_iter=dav_iter, &
                             dgl_shift=shift, &
                             dgl_tol=tol, &
                             dgl_memory=memory, &
                             dgl_memory_unit=memory_unit &
                             )
!
    if (ok) then
        write (6, *) ' Davidson converged.'
        call dump_eigpairs(lutest, n, n_targ, eig, evec, "Davidson")
    else
        write (6, *) ' Davidson failed to converge.'
    end if
!
! test generalized davidson:
!
    call init_eigenpairs(n, n_max, eig, evec, ok)
!
    write (6, *) ' testing Generalized Davidson:'
    call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
                             dgl_verbose=verbose, &
                             dgl_max_iter=max_iter, &
                             dgl_dav_iter=dav_iter, &
                             dgl_shift=shift, &
                             dgl_tol=tol, &
                             dgl_memory=memory, &
                             dgl_memory_unit=memory_unit &
                             )
!
    if (ok) then
        write (6, *) ' Generalized Davidson converged.'
        call dump_eigpairs(lutest, n, n_targ, eig, evec, "Generalized Davidson")
    else
        write (6, *) 'Generalized Davidson failed to converge.'
    end if
!
! test non-symmetric davidson:
!
    call init_eigenpairs(n, n_max, eig, evec, ok)
    call init_eigenpairs(n, n_max, eig, evec_2, ok)
!
    write (6, *) ' testing non-symmetric Davidson:'
    call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "LR", eig, evec, ok, evec_2=evec_2, &
                                   dgl_verbose=verbose, &
                                   dgl_max_iter=max_iter, &
                                   dgl_dav_iter=dav_iter, &
                                   dgl_shift=shift, &
                                   dgl_tol=tol, &
                                   dgl_memory=memory, &
                                   dgl_memory_unit=memory_unit &
                                   )
!
    if (ok) then
        write (6, *) ' non-symmetric Davidson converged.'
        call dump_eigpairs(lutest, n, n_targ, eig, evec, "Non Symmetric Davidson, Right")
        call dump_eigpairs(lutest, n, n_targ, eig, evec_2, "Non Symmetric Davidson, Left")
    else
        write (6, *) ' non-symmetric Davidson failed to converge.'
    end if
!
! test lobpcg:
!
    call init_eigenpairs(n, n_max, eig, evec, ok)
!
    write (6, *) ' testing LOBPCG:'

    call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                           dgl_verbose=verbose, &
                           dgl_max_iter=max_iter, &
                           dgl_shift=shift, &
                           dgl_tol=tol, &
                           dgl_memory=memory, &
                           dgl_memory_unit=memory_unit &
                           )
!
    if (ok) then
        write (6, *) ' LOBPCG converged.'
        call dump_eigpairs(lutest, n, n_targ, eig, evec, "LOBPCG")
    else
        write (6, *) ' LOBPCG failed to converge.'
    end if
!
! test generalized lobpcg:
!
    call init_eigenpairs(n, n_max, eig, evec, ok)
!
    write (6, *) ' testing Generalized LOBPCG:'

    call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
                           dgl_verbose=verbose, &
                           dgl_max_iter=max_iter, &
                           dgl_shift=shift, &
                           dgl_tol=tol, &
                           dgl_memory=memory, &
                           dgl_memory_unit=memory_unit &
                           )
!
    if (ok) then
        write (6, *) ' Generalized LOBPCG converged.'
        call dump_eigpairs(lutest, n, n_targ, eig, evec, "Generalized LOBPCG")
    else
        write (6, *) ' Generalized LOBPCG failed to converge.'
    end if
!
! test smogd:
!
    deallocate (evec)
    allocate (evec(2*n, n_max))
    call init_eigenpairs(2*n, n_max, eig, evec, ok)
!
    write (6, *) ' testing SMOGD:'
    call dgl_smogd_driver(2*n, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, &
                          eig, evec, ok, &
                          dgl_verbose=verbose, &
                          dgl_max_iter=max_iter, &
                          dgl_dav_iter=dav_iter, &
                          dgl_tol=tol, &
                          dgl_memory=memory, &
                          dgl_memory_unit=memory_unit &
                          )
    if (ok) then
        write (6, *) ' SMOGD converged.'
        call dump_eigpairs(lutest, n*2, n_targ, eig, evec, "SMOGD")
    else
        write (6, *) ' SMOGD failed to converge.'
    end if
!
! close the output file:
!
    close (lutest)
!
! free the memory:
!
    deallocate (evec, evec_2, eig)
!
end program test_fortran

