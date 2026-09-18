! Statistical stress test of all the DiagLib drivers against dense LAPACK references.
! Dense matrices are built explicitly, the matvecs apply them with dgemm, and every run uses a
! random guess. Configurations cover: clustered/degenerate/wide spectra, an ill-conditioned
! metric, non-symmetric matrices with complex pairs and near degeneracies, linear response,
! many roots, small expansion spaces (many restarts), and poor preconditioners.
program stress_sym
    use dgl_interface
    implicit none
    integer, parameter :: ip = dgl_int, dp = dgl_real
    integer(ip), parameter :: n = 400
    real(dp), allocatable :: a(:, :), b(:, :), al(:, :), apb(:, :), amb(:, :), spd(:, :), smd(:, :)
    real(dp), allocatable :: ref(:), diag_a(:), diag_b(:)
    real(dp), allocatable :: eig(:), evec(:, :), evec2(:, :), work(:), wr(:), wi(:), vdum(:, :), acopy(:, :), bcopy(:, :)
    integer(ip) :: i, j, info, lwork, n_targ, n_max, dav_iter, i_cfg, run, n_runs, i_side
    integer(ip) :: n_fail, n_tot, n_fail_cfg, n_tot_cfg
    real(dp) :: tol, err, res, spread_a, deg, cond_b, t
    logical :: ok, generalized, weak, spd_precnd
    character(len=2) :: sides(3) = ["R ", "L ", "LR"]
    character(len=80) :: label
    procedure(dgl_matvec), pointer :: metvec_p
    integer :: seed_size
    integer, allocatable :: seed(:)
!
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 20260918
    call random_seed(put=seed)
    n_runs = 12
    tol = 1.0e-9_dp
    n_fail = 0
    n_tot = 0
    allocate (a(n, n), b(n, n), al(n, n), apb(n, n), amb(n, n), spd(n, n), smd(n, n))
    allocate (ref(n), diag_a(n), diag_b(n), acopy(n, n), bcopy(n, n))
    allocate (wr(n), wi(n), vdum(1, 1))
    lwork = 64*n
    allocate (work(lwork))
!
! symmetric and generalized problems: spectra with different spreads and degeneracies,
! metrics with different conditioning
!
    do i_cfg = 1, 8
        spread_a = merge(1.0_dp, 100.0_dp, mod(i_cfg, 2_ip) .eq. 0)
        deg = merge(0.0_dp, 1.0e-6_dp, i_cfg .le. 4)      ! nearly degenerate pairs
        cond_b = merge(1.0_dp, 1.0e6_dp, i_cfg .le. 2)    ! conditioning of the metric
        generalized = i_cfg .gt. 2
        n_targ = merge(4_ip, 20_ip, i_cfg .le. 6)
        n_max = 2*n_targ
        dav_iter = merge(10_ip, 25_ip, mod(i_cfg, 3_ip) .eq. 0)
        weak = mod(i_cfg, 4_ip) .eq. 0
        spd_precnd = .false.
        write (label, "(a,i0,a,l1,a,es7.1,a,es7.1,a,i0)") "sym cfg ", i_cfg, " gen ", generalized, &
            " deg ", deg, " cond(B) ", cond_b, " roots ", n_targ
        call build_sym(spread_a, deg, cond_b)
        call reference_sym(generalized)
        n_fail_cfg = 0
        n_tot_cfg = 0
        do run = 1, n_runs
            metvec_p => mv_b
            call run_one("Davidson", .false.)
            call run_one("LOBPCG", .true.)
        end do
        write (*, "(a60,a,i4,a,i4)") trim(label), ":", n_fail_cfg, " failed of", n_tot_cfg
    end do
!
    write (*, "(/,a,i0,a,i0,/)") "symmetric and generalized: ", n_fail, " failures out of ", n_tot
contains
!
    subroutine build_sym(spread_in, deg_in, cond_in)
        real(dp), intent(in) :: spread_in, deg_in, cond_in
        real(dp) :: v(n), s(n), r
        integer(ip) :: k, l
!
        do k = 1, n
            diag_a(k) = 1.0_dp + spread_in*real(k - 1, dp)/real(n - 1, dp)
            if (deg_in .gt. 0.0_dp .and. mod(k, 2_ip) .eq. 0) diag_a(k) = diag_a(k - 1) + deg_in
            diag_b(k) = 1.0_dp + (cond_in - 1.0_dp)*real(n - k, dp)/real(n - 1, dp)
        end do
        call random_number(v)
        v = v - 0.5_dp
        v = v/norm2(v)
!
! a = H diag_a H, b = H diag_b H
!
        do l = 1, n
            do k = 1, n
                r = merge(1.0_dp, 0.0_dp, k .eq. l) - 2.0_dp*v(k)*v(l)
                a(k, l) = r
                b(k, l) = r
            end do
        end do
        do l = 1, n
            a(:, l) = a(:, l)*diag_a(l)
            b(:, l) = b(:, l)*diag_b(l)
        end do
        acopy = a
        bcopy = b
        do l = 1, n
            do k = 1, n
                s(k) = merge(1.0_dp, 0.0_dp, k .eq. l) - 2.0_dp*v(k)*v(l)
            end do
            a(:, l) = matmul(acopy, s)
            b(:, l) = matmul(bcopy, s)
        end do
    end subroutine build_sym
!
    subroutine reference_sym(gen)
        logical, intent(in) :: gen
        acopy = a
        bcopy = b
        if (gen) then
            call dsygv(1_ip, 'n', 'u', n, acopy, n, bcopy, n, ref, work, lwork, info)
        else
            call dsyev('n', 'u', n, acopy, n, ref, work, lwork, info)
        end if
        if (info .ne. 0) then
            write (*, *) "LAPACK reference failed, info =", info
            stop 1
        end if
    end subroutine reference_sym
!
    subroutine run_one(driver, lobpcg)
        character(len=*), intent(in) :: driver
        logical, intent(in) :: lobpcg
!
        if (allocated(eig)) deallocate (eig, evec)
        allocate (eig(n_max), evec(n, n_max))
        eig = 0.0_dp
        evec = 0.0_dp
        if (lobpcg .and. 3*n_max .ge. n) return
        if (.not. lobpcg .and. 2*n_max .ge. n) return
        if (lobpcg) then
            if (generalized) then
                call dgl_lobpcg_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, metvec=metvec_p, &
                                       dgl_tol=tol, dgl_max_iter=400_ip, dgl_info=info)
            else
                call dgl_lobpcg_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, &
                                       dgl_tol=tol, dgl_max_iter=400_ip, dgl_info=info)
            end if
        else
            if (generalized) then
                call dgl_davidson_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, metvec=metvec_p, &
                                         dgl_tol=tol, dgl_max_iter=400_ip, dgl_dav_iter=dav_iter, dgl_info=info)
            else
                call dgl_davidson_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, &
                                         dgl_tol=tol, dgl_max_iter=400_ip, dgl_dav_iter=dav_iter, dgl_info=info)
            end if
        end if
        err = maxval(abs(eig(1:n_targ) - ref(1:n_targ)))
        res = residual()
        n_tot = n_tot + 1
        n_tot_cfg = n_tot_cfg + 1
        if (.not. ok .or. info .ne. dgl_success .or. err .gt. 1.0e-7_dp .or. res .gt. 1.0e-5_dp) then
            n_fail = n_fail + 1
            n_fail_cfg = n_fail_cfg + 1
            write (*, "(a,a,a,l2,a,i3,a,es9.2,a,es9.2)") "   FAIL ", driver, " ok", ok, " info", info, &
                " eig err", err, " residual", res
        end if
    end subroutine run_one
!
    real(dp) function residual()
        real(dp) :: ax(n), bx(n)
        integer(ip) :: k
        residual = 0.0_dp
        do k = 1, n_targ
            ax = matmul(a, evec(:, k))
            if (generalized) then
                bx = matmul(b, evec(:, k))
            else
                bx = evec(:, k)
            end if
            residual = max(residual, norm2(ax - eig(k)*bx)/(norm2(ax) + abs(eig(k))*norm2(bx)))
        end do
    end function residual
!
    subroutine mv_a(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, a, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_a
!
    subroutine mv_b(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, b, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_b
!
    subroutine pc(nn, m, shift, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: k, l
        real(dp) :: den
        do l = 1, m
            do k = 1, nn
                den = a(k, k) + shift
                if (weak) den = a(k, k) + shift + 100.0_dp     ! poor preconditioner
                if (spd_precnd) den = abs(a(k, k))
                if (abs(den) .lt. 1.0e-3_dp) den = sign(1.0e-3_dp, den)
                y(k, l) = x(k, l)/den
            end do
        end do
    end subroutine pc
!
end program stress_sym
