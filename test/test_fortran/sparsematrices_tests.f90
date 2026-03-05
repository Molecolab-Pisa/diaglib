module sparsematrices
use utility
implicit none

private

real(dp), allocatable, dimension(:,:) :: matrix
!! container for the matrix
integer :: nrows
!! Rows number
integer :: ncols
!! Columns number
integer :: nonzero
!! Number of non-zero elements
integer, parameter :: useless_lines = 13
!! Useless lines in matrix files

integer :: i, j, k

character(len=30), public :: matrix_file(2) = [ "494_bus.mtx ", &
                                        "1138_bus.mtx" ]
!! Files containing the matrices to diagonalize

public read_dimensions, read_matrix

contains

subroutine read_dimensions(fname)
!! Read dimensions for a matrix
implicit none
character(len=*) :: fname
integer :: unit

open(file=fname, newunit=unit, status="old")

do i = 1, useless_lines
    read(unit,*)
enddo

read(unit,*) nrows, ncols, nonzero
close(unit)

end subroutine read_dimensions

subroutine read_matrix(fname)
!! Allocate and Read a matrix in a global record that can be used for a matvec
implicit none
character(len=*) :: fname
integer :: unit
real(dp) :: val

allocate(matrix(nrows,ncols))

open(file=fname, newunit=unit, status="old")

do i = 1, useless_lines + 1
    read(unit,*)
enddo

do k = 1, nonzero
    read(unit,*) i, j, val
    matrix(i,j) = val
enddo
close(unit)

end subroutine read_matrix

subroutine matvec()
!! Perform a matvec
end subroutine matvec

subroutine compare_exact_sym()
!! Diagonalize exactly the simmetric matrix and compare with the solver results
end subroutine compare_exact_sym

subroutine compare_exact_nosym()
!! Diagonalize exactly the non-symmetric matrix and compare with the solver results
end subroutine compare_exact_nosym



end module sparsematrices