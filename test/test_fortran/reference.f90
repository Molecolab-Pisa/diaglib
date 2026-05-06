program reference
!!
!! Build reference file for the results of DiagLib tests
!!
    use utility
    implicit none

    integer :: lwork = 10000, ilwork = 10000, info
    real(dp), allocatable :: a(:, :), b(:, :), copy(:, :), cont(:, :)
    real(dp), allocatable :: eig(:, :)
    real(dp), allocatable :: evec(:, :, :)
    real(dp), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    allocate (a(n, n), b(n, n), copy(n, n), cont(n, n))
    allocate (eig(n, 2), evec(n, n, 2))
    allocate (work(lwork), iwork(ilwork))
    a = zero; b = zero

    open (file=trim(reference_fname), unit=luref, status="unknown")

    write (*, f_string) "Building symmetric matrix"
    call get_sym_matrix(n, a)
    write (*, f_string) "Building symmetric metric"
    call get_sym_metric(n, b)
    copy = a

    !symmetric diagonalization
    write (*, f_string) "Running symmetric diagonalization"
    call dsyev("V", "U", n, a, n, eig, work, lwork, info)
    call check_lapack(info)
    call dump_eigpairs(luref, n, n, eig, a, "Symmetric diagonalization")

    a = copy
    !symmetric generalized problem
    write (*, f_string) "Running symmetric generalized diagonalization"
    call dsygv(1, "V", "U", n, a, n, b, n, eig, work, lwork, info)
    call check_lapack(info)
    call dump_eigpairs(luref, n, n, eig, a, "Symmetric Generalized diagonalization")

    write (*, f_string) "Building non symmetric matrix"
    call get_asym_matrix(n, a)

    !non symmetric diagonalization
    write (*, f_string) "Running non symmetric generalized diagonalization"
    call dgeev("V", "V", n, a, n, eig, eig(1, 2), evec, n, evec(1, 1, 2), n, work, lwork, info)
    call check_lapack(info)

    if (sqrt(dot_product(eig(:, 2), eig(:, 2))) .gt. 1.e-12_dp) write (*, f_string) "Immaginary eigs detected !!"
    call sort_eigenpairs(n, n, eig, evec)

    call dump_eigpairs(luref, n, n, eig, evec, "Non Symmetric diagonalization, Right")
    call dump_eigpairs(luref, n, n, eig, evec(1, 1, 2), "Non Symmetric diagonalization, Left")

    deallocate (a, b)
    deallocate (eig)

    allocate (a(2*n, 2*n), b(2*n, 2*n))
    allocate (eig(2*n, 1))
    a = zero; b = zero

    write (*, f_string) "Building linear response matrices"
    call get_apb_matrix(n, copy)
    call get_amb_matrix(n, cont)

    a(1:n, 1:n) = (copy + cont)*half
    a(n + 1:2*n, n + 1:2*n) = (copy + cont)*half
    a(1:n, n + 1:2*n) = (copy - cont)*half
    a(n + 1:2*n, 1:n) = (copy - cont)*half

    call get_spd_matrix(n, copy)
    call get_smd_matrix(n, cont)

    b(1:n, 1:n) = (copy + cont)*half
    b(n + 1:2*n, n + 1:2*n) = -(copy + cont)*half
    b(1:n, n + 1:2*n) = (copy - cont)*half
    b(n + 1:2*n, 1:n) = -(copy - cont)*half

    !linear response problem
    write (*, f_string) "Running linear response diagonalization"
    call dsygv(1, "V", "L", 2*n, b, 2*n, a, 2*n, eig, work, lwork, info)
    call check_lapack(info)

    eig = -one/eig
    do i = 1, 2*n
        a(:, i) = b(:, 2*n + 1 - i)
    end do
    call dump_eigpairs(luref, 2*n, 2*n, eig, a, "Linear response diagonalization")

    write (*, f_string) "All done!"
    deallocate (evec, eig)
    deallocate (a, b, copy, cont)
    deallocate (work)

    close (luref)

end program reference

subroutine get_sym_matrix(n, mat)

    use utility, only: dp, i, j, one
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, i - 1
            mat(j, i) = one/real(i + j, dp)
            !mat(i,j) = mat(j,i)
        end do
        mat(i, i) = one + real(i, dp)
    end do

end subroutine get_sym_matrix

subroutine get_asym_matrix(n, mat)

    use utility, only: dp, i, j, one
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, n
            mat(j, i) = real(i, dp)/real(j*(i + j), dp)
        end do
        mat(i, i) = one + real(i, dp)
    end do

end subroutine get_asym_matrix

subroutine get_sym_metric(n, mat)

    use utility, only: dp, i, j, one
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, i - 1
            mat(j, i) = one/real(i + j, dp)
            !mat(i,j) = mat(j,i)
        end do
        mat(i, i) = one
    end do

end subroutine get_sym_metric

subroutine get_apb_matrix(n, mat)

    use utility, only: dp, i, j, one, five
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, i - 1
            mat(j, i) = one/real(i + j, dp)
            mat(i, j) = mat(j, i)
        end do
        mat(i, i) = five + real(i, dp)
    end do

end subroutine get_apb_matrix

subroutine get_amb_matrix(n, mat)

    use utility, only: dp, i, j, two
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, i - 1
            mat(j, i) = 0.20_dp/real(i + j, dp)
            mat(i, j) = mat(j, i)
        end do
        mat(i, i) = two + real(i, dp)
    end do

end subroutine get_amb_matrix

subroutine get_spd_matrix(n, mat)

    use utility, only: dp, i, j, one
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, i - 1
            mat(j, i) = -0.05_dp
        end do
        do j = i + 1, n
            mat(j, i) = +0.05_dp
        end do
        mat(i, i) = one
    end do

end subroutine get_spd_matrix

subroutine get_smd_matrix(n, mat)

    use utility, only: dp, i, j, one
    implicit none
    integer, intent(in) :: n
    real(dp) :: mat(n, n)
    do i = 1, n
        do j = 1, i - 1
            mat(j, i) = +0.05_dp
        end do
        do j = i + 1, n
            mat(j, i) = -0.05_dp
        end do
        mat(i, i) = one
    end do

end subroutine get_smd_matrix

subroutine check_lapack(info)
    implicit none
    integer, intent(in) :: info

    if (info .ne. 0) then
        write (*, "(t3,a)") "Lapack Failed"
        stop 1
    end if

end subroutine check_lapack
