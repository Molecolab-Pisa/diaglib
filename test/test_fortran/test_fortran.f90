program test_fortran
!!
!! tests the various functionalities of diaglib.
!!
    use dgl_interface
    use direct_matvecs
    use utility
    implicit none
    logical :: ok
    real(dp), allocatable :: eig(:), evec(:, :), evec_2(:, :)
!
! Check if the reference exists
!
    call check_reference(n)
!
! open a text file to dump the output.
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
!====================================
! TESTING AREA
!====================================
!
    call test_davidson
    call test_davidson_generalized
    call test_lobpcg
    call test_lobpcg_generalized
    call test_nonsym_davidson
    call test_smogd
!
! close the output file:
!
    close (lutest)
!
! free the memory:
!
    deallocate (evec, evec_2, eig)
!    
contains
!
    subroutine test_davidson()
        implicit none
!
! test davidson:
!
        call init_eigenpairs(n, n_max, eig, evec, ok)
!
        write (6, f_string) 'testing Davidson:'
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
            write (6, f_string) 'Davidson converged.'
            eig = eig - shift
            ok = compare_eigs(n, n_targ, tol, eig, evec, "Symmetric diagonalization")
            call dump_eigpairs(lutest, n, n_targ, eig, evec, "Davidson")
            if (.not. ok) then
                write (*, f_string) "Davidson results do not match with reference! Maybe rerun reference?"
            end if
        else
            write (6, f_string) 'Davidson failed to converge.'
        end if
!
    end subroutine test_davidson
!
    subroutine test_davidson_generalized()
        implicit none
!
! test generalized davidson:
!
        call init_eigenpairs(n, n_max, eig, evec, ok)
!
        write (6, f_string) 'testing Generalized Davidson:'
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
            write (6, f_string) 'Generalized Davidson converged.'
            ok = compare_eigs(n, n_targ, tol, eig, evec, "Symmetric Generalized")
            call dump_eigpairs(lutest, n, n_targ, eig, evec, "Generalized Davidson")
            if (.not. ok) then
                write (*, f_string) "Generalized Davidson results do not match with reference! Maybe rerun reference?"
            end if

        else
            write (6, f_string) 'Generalized Davidson failed to converge.'
        end if
!
    end subroutine test_davidson_generalized
!
    subroutine test_nonsym_davidson()
        implicit none
!
! test non-symmetric davidson:
!
        call init_eigenpairs(n, n_max, eig, evec, ok)
        call init_eigenpairs(n, n_max, eig, evec_2, ok)
!
        write (6, f_string) 'testing non-symmetric Davidson:'
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
            write (6, f_string) 'non-symmetric Davidson converged.'
            eig = eig - shift
            ok = compare_eigs(n, n_targ, tol, eig, evec, "Non Symmetric diagonalization, Right")
            call dump_eigpairs(lutest, n, n_targ, eig, evec, "Non Symmetric Davidson, Right")
            if (.not. ok) then
                write (*, f_string) "Right Non symmetric Davidson results do not match with reference! Maybe rerun reference?"
            end if

            ok = compare_eigs(n, n_targ, tol, eig, evec_2, "Non Symmetric diagonalization, Left")
            call dump_eigpairs(lutest, n, n_targ, eig, evec_2, "Non Symmetric Davidson, Left")
            if (.not. ok) then
                write (*, f_string) "Left Non symmetric Davidson results do not match with reference! Maybe rerun reference?"
            end if
        else
            write (6, f_string) 'non-symmetric Davidson failed to converge.'
        end if
!
    end subroutine test_nonsym_davidson
!
    subroutine test_lobpcg()
        implicit none

!
! test lobpcg:
!
        call init_eigenpairs(n, n_max, eig, evec, ok)
!
        write (6, f_string) 'testing LOBPCG:'

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
            write (6, f_string) 'LOBPCG converged.'
            eig = eig - shift
            ok = compare_eigs(n, n_targ, tol, eig, evec, "Symmetric diagonalization")
            call dump_eigpairs(lutest, n, n_targ, eig, evec, "LOBPCG")
            if (.not. ok) then
                write (*, f_string) "LOBPCG results do not match with reference! Maybe rerun reference?"
            end if
        else
            write (6, f_string) 'LOBPCG failed to converge.'
        end if

    end subroutine test_lobpcg
!
    subroutine test_lobpcg_generalized()
        implicit none
!
! test generalized lobpcg:
!
        call init_eigenpairs(n, n_max, eig, evec, ok)
!
        write (6, f_string) 'testing Generalized LOBPCG:'

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
            write (6, f_string) 'Generalized LOBPCG converged.'
            ok = compare_eigs(n, n_targ, tol, eig, evec, "Symmetric Generalized")
            call dump_eigpairs(lutest, n, n_targ, eig, evec, "Generalized LOBPCG")
            if (.not. ok) then
                write (*, f_string) "Generalized LOBPCG results do not match with reference! Maybe rerun reference?"
            end if
        else
            write (6, f_string) 'Generalized LOBPCG failed to converge.'
        end if
!
    end subroutine test_lobpcg_generalized
!
    subroutine test_smogd()
        implicit none
!
! test smogd:
!
        deallocate (evec, eig)
        allocate (evec(2*n, n_max), eig(2*n))
        call init_eigenpairs(2*n, n_max, eig, evec, ok)
!
        write (6, f_string) 'testing SMOGD:'
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
            write (6, f_string) 'SMOGD converged.'
            ok = compare_eigs(2*n, n_targ, tol, eig, evec, "Linear response")
            call dump_eigpairs(lutest, n*2, n_targ, eig, evec, "SMOGD")
            if (.not. ok) then
                write (*, f_string) "SMOGD results do not match with reference! Maybe rerun reference?"
            end if
        else
            write (6, f_string) 'SMOGD failed to converge.'
        end if
!
        deallocate (evec, eig)
        allocate (evec(n, n_max), eig(n))
!
    end subroutine test_smogd
!
end program test_fortran

