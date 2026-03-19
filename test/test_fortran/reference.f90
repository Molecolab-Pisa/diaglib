program reference
!!
!! Build reference file for the results of DiagLib tests
!!
    use utility
    implicit none

    integer :: lwork = 10000, ilwork = 10000, info
    real(dp), allocatable :: a(:, :), b(:, :), copy(:, :)
    real(dp), allocatable :: eig(:, :)
    real(dp), allocatable :: evec(:, :, :)
    real(dp), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    allocate (a(n, n), b(n, n), copy(n, n))
    allocate (eig(n, 2), evec(n, n, 2))
    allocate (work(lwork), iwork(ilwork))
    a = zero; b = zero

    open (file=trim(reference_fname), unit=luref, status="unknown")

    call get_sym_matrix(n, a)
    call get_sym_metric(n, b)
    copy = a

    !symmetric diagonalization
    call dsyev("V", "U", n, a, n, eig, work, lwork, info)
    call check_lapack(info)
    call dump_eigpairs(luref, n, n_targ, eig, a, "Symmetric diagonalization")

    a = copy
    !symmetric generalized problem
    call dsygv(1, "V", "U", n, a, n, b, n, eig, work, lwork, info)
    call check_lapack(info)
    call dump_eigpairs(luref, n, n_targ, eig, a, "Symmetric Generalized diagonalization")

    call get_asym_matrix(n, a)
    !non symmetric diagonalization
    call dgeev("V", "V", n, a, n, eig, eig(1, 2), evec, n, evec(1, 1, 2), n, work, lwork, info)
    call check_lapack(info)

    if (sqrt(dot_product(eig(:, 2), eig(:, 2))) .gt. 1.e-12_dp) print *, "Immaginary eigs detected"
    call sort_eigenpairs(n, n_targ, eig, evec)

    call dump_eigpairs(luref, n, n_targ, eig, evec, "Non Symmetric diagonalization, Left")
    call dump_eigpairs(luref, n, n_targ, eig, evec(1, 1, 2), "Non Symmetric diagonalization, Right")

    deallocate (a, b)
    deallocate (eig)

    allocate (a(2*n, 2*n), b(2*n, 2*n))
    allocate (eig(2*n, 1))
    a = zero; b = zero

    call get_apb_matrix(n, copy)
    a(1:n, 1:n) = copy

    call get_amb_matrix(n, copy)
    a(n + 1:2*n, n + 1:2*n) = copy

    call get_spd_matrix(n, copy)
    !b(n + 1:2*n, 1:n) = -copy
    b(1:n, n + 1:2*n) = copy
    
    call get_smd_matrix(n, copy)
    !b(1:n, n + 1:2*n) = -copy
    b(n + 1:2*n, 1:n) = copy

    call prtmat(2*n, 2*n, a)
    print *
    call prtmat(2*n, 2*n, b)

    !linear response problem
    call dsygv(1, "V", "U", 2*n, b, 2*n, a, 2*n, eig, work, lwork, info)
    call check_lapack(info)
    eig = -one/eig
    call dump_eigpairs(luref, 2*n, n_targ, eig, evec, "Linear response diagonalization")

    deallocate (evec, eig)
    deallocate (a, b, copy)
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
            mat(j, i) = +0.05_dp
        end do
        do j = i + 1, n
            mat(j, i) = -0.05_dp
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
            mat(j, i) = -0.05_dp
        end do
        do j = i + 1, n
            mat(j, i) = +0.05_dp
        end do
        mat(i, i) = one
    end do

end subroutine get_smd_matrix

subroutine check_lapack(info)
    implicit none
    integer, intent(in) :: info

    if (info .ne. 0) then
        stop "Lapack Failed"
    end if

end subroutine check_lapack
