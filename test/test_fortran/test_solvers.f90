module solvers
    use dgl_interface
    use matvecs

    integer :: i, j
    logical :: ok
    integer, parameter :: lutest = 100
    real(dgl_real), allocatable :: eig(:), evec(:, :), evec_2(:, :)
    integer, parameter :: memory = 100
    character(len=2), parameter :: memory_unit = "MB"
    procedure(), pointer :: mx_p => null()

contains

    subroutine reset_output_file
        implicit none

        open(file="output_fortran.txt",unit=lutest,status="unknown",position="rewind")
        close(lutest)
        
    end subroutine reset_output_file

    subroutine test_davidson(n, n_targ, n_max, verbose, max_iter, dav_iter, tol, shift)
!!
!! Test davidson:
!!
        implicit none
        integer, intent(in) :: n, n_targ, n_max
        integer, intent(in), optional :: max_iter, dav_iter
        real(dgl_real), intent(in), optional :: tol, shift
        logical, intent(in), optional :: verbose
!
        call init_eigenpairs(n, n_max)
!
        write (6, *) ' testing Davidson:'
        call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                                 dgl_verbose=verbose, &
                                 dgl_max_iter=max_iter, &
                                 dgl_dav_iter=dav_iter, &
                                 dgl_tol=tol, &
                                 dgl_shift=shift, &
                                 dgl_memory=memory, &
                                 dgl_memory_unit=memory_unit)
!
        if (ok) then
            write (6, *) ' Davidson converged.'
            call dump_eigpairs(n, n_targ)
        else
            write (6, *) ' Davidson failed to converge.'
        end if
!
        deallocate (eig, evec)
!
    end subroutine test_davidson
!
    subroutine test_lobpcg(n, n_targ, n_max, verbose, max_iter, tol)
!!
!! Test lobpcg:
!!
        implicit none
        integer, intent(in) :: n, n_targ, n_max
        integer, intent(in), optional :: max_iter
        real(dgl_real), intent(in), optional :: tol
        logical, intent(in), optional :: verbose
!
        call init_eigenpairs(n, n_max)
!
        write (6, *) ' testing LOBPCG:'
        call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
                                 dgl_verbose=verbose, &
                                 dgl_max_iter=max_iter, &
                                 dgl_tol=tol, &
                                 dgl_memory=memory, &
                                 dgl_memory_unit=memory_unit)
!
        if (ok) then
            write (6, *) ' LOBPCG converged.'
            call dump_eigpairs(n, n_targ)
        else
            write (6, *) ' LOBPCG failed to converge.'
        end if
!
        deallocate (eig, evec)
!
    end subroutine test_lobpcg
!    
    subroutine test_nosym_davidson(n, n_targ, n_max, side, verbose, max_iter, dav_iter, tol, shift)
!!
!! Test davidson:
!!
        implicit none
        integer, intent(in) :: n, n_targ, n_max
        integer, intent(in), optional :: max_iter, dav_iter
        character(len=2), intent(in) :: side
        real(dgl_real), intent(in), optional :: tol, shift
        logical, intent(in), optional :: verbose
!
        call init_eigenpairs(n, n_max)
!
        write (6, *) ' testing Davidson:'
        call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, side, eig, evec, ok, evec_2 = evec_2, &
                                 dgl_verbose=verbose, &
                                 dgl_max_iter=max_iter, &
                                 dgl_dav_iter=dav_iter, &
                                 dgl_tol=tol, &
                                 dgl_shift=shift, &
                                 dgl_memory=memory, &
                                 dgl_memory_unit=memory_unit)
!
        if (ok) then
            write (6, *) ' Davidson converged.'
            call dump_eigpairs(n, n_targ)
        else
            write (6, *) ' Davidson failed to converge.'
        end if
!
        deallocate (eig, evec)
!
    end subroutine test_nosym_davidson
!
    subroutine test_smogd(n, n_targ, n_max, verbose, max_iter, dav_iter, tol)
!!
!! Test smogd:
!!
        implicit none
        integer, intent(in) :: n, n_targ, n_max
        integer, intent(in), optional :: max_iter, dav_iter
        real(dgl_real), intent(in), optional :: tol
        logical, intent(in), optional :: verbose
!
        call init_eigenpairs(n, n_max)
!
        write (6, *) ' testing SMOGD:'
        call dgl_smogd_driver(2*n, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, eig, evec, ok, &
                                 dgl_verbose=verbose, &
                                 dgl_max_iter=max_iter, &
                                 dgl_dav_iter=dav_iter, &
                                 dgl_tol=tol, &
                                 dgl_memory=memory, &
                                 dgl_memory_unit=memory_unit)
!
        if (ok) then
            write (6, *) ' SMOGD converged.'
            call dump_eigpairs(n, n_targ)
        else
            write (6, *) ' SMOGD failed to converge.'
        end if
!
        deallocate (eig, evec)
!
    end subroutine test_smogd
!
    subroutine init_eigenpairs(n, n_max)
!!
!! Make a simple guess
!!   
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: n_max
!
        allocate (eig(n_max), evec(n, n_max))
        ok = .false.
        eig = dgl_zero
        evec = dgl_zero
        do i = 1, n_max
            evec(i, i) = dgl_one
        end do
    end subroutine init_eigenpairs

    subroutine dump_eigpairs(n, n_targ)
!!
!! open a text file for the output, to be used to compare the results with a reference.
!!
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: n_targ
        integer :: i, j
        open (file='output_fortran.txt', form='formatted', &
              unit=lutest, status='old', position='append')
        !write (lutest, 1000) 'Davidson'
        write (lutest, *)
        write (lutest, 1010)
        do i = 1, n_targ
            write (lutest, 1020) i, eig(i)
        end do
        do i = 1, n_targ
            if (evec(1, i) .lt. dgl_zero) evec(:, i) = -evec(:, i)
        end do
        write (lutest, 1031) (i, i=1, n_targ)
        do j = 1, n
            write (lutest, 1022) j, (evec(j, i), i=1, n_targ)
        end do
        write (lutest, *)
        close (lutest)
!
1000    format(t3, a, 1x, 'results')
1010    format(t3, 'Eigenvalues:')
1020    format(t3, i5, f14.6)
1021    format(t3, i5, f12.4)
1022    format(t3, i5, *(f14.6))
1030    format(t3, 'Eigenvector ', i3, ':')
1031    format(t3, 'Eigenvector ', tl1, *(i8, 6x))
1032    format(t3, a6, 'Eigenvector ', tl7, *(i8, 6x))
!
    end subroutine dump_eigpairs
!
end module solvers
