! Stress test of the non-symmetric and linear-response drivers against dense LAPACK references.
! Non-symmetric: M = P diag P^-1 with P = I + u w^T, spectra with nearly degenerate pairs and
! complex pairs among the lowest roots. Linear response: A+B, A-B positive definite, S+D, S-D
! with an antisymmetric part. Every run uses a random guess; small expansion spaces force
! restarts.
program stress_nosym
    use dgl_interface
    implicit none
    integer, parameter :: ip = dgl_int, dp = dgl_real
    integer(ip), parameter :: n = 300
    real(dp) :: m_r(n, n), m_l(n, n), apb(n, n), amb(n, n), spd(n, n), smd(n, n)
    real(dp) :: e2(2*n, 2*n), met(2*n, 2*n), acopy(2*n, 2*n), bcopy(2*n, 2*n)
    real(dp) :: ref(n), ref_lr(2*n), wr(2*n), wi(2*n), vdum(1, 1)
    real(dp), allocatable :: work(:), eig(:), evec(:, :), evec2(:, :), evlr(:, :)
    integer(ip) :: i, j, k, info, lwork, n_targ, n_max, dav_iter, i_cfg, run, i_side, n_kept
    integer(ip) :: n_fail, n_tot, n_fail_cfg, n_tot_cfg
    real(dp) :: tol, err, res, deg, cplx
    logical :: ok
    character(len=2) :: sides(3) = ["R ", "L ", "LR"]
    character(len=90) :: label
    integer :: seed_size
    integer, allocatable :: seed(:)
!
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 777
    call random_seed(put=seed)
    lwork = 64*2*n
    allocate (work(lwork))
    tol = 1.0e-9_dp
    n_fail = 0
    n_tot = 0
!
! non-symmetric driver
!
    do i_cfg = 1, 6
        deg = merge(0.0_dp, 1.0e-5_dp, i_cfg .le. 2)
        cplx = merge(0.0_dp, 0.3_dp, i_cfg .ge. 5)
        n_targ = merge(4_ip, 12_ip, i_cfg .le. 4)
        n_max = 2*n_targ
        dav_iter = merge(10_ip, 25_ip, mod(i_cfg, 2_ip) .eq. 0)
        call build_nosym(deg, cplx)
        call reference_nosym()
        n_fail_cfg = 0
        n_tot_cfg = 0
        write (label, "(a,i0,a,es8.1,a,es8.1,a,i0,a,i0)") "nosym cfg ", i_cfg, " deg ", deg, " cplx ", cplx, &
            " roots ", n_targ, " dav_iter ", dav_iter
        do run = 1, 6
            do i_side = 1, 3
                call run_nosym(sides(i_side))
            end do
        end do
        write (*, "(a70,a,i4,a,i4)") trim(label), ":", n_fail_cfg, " failed of", n_tot_cfg
    end do
!
! linear response driver
!
    do i_cfg = 1, 4
        n_targ = merge(4_ip, 10_ip, i_cfg .le. 2)
        n_max = 2*n_targ
        dav_iter = merge(10_ip, 25_ip, mod(i_cfg, 2_ip) .eq. 0)
        call build_lr(merge(0.0_dp, 0.2_dp, i_cfg .le. 2))
        call reference_lr()
        n_fail_cfg = 0
        n_tot_cfg = 0
        write (label, "(a,i0,a,i0,a,i0)") "smogd cfg ", i_cfg, " roots ", n_targ, " dav_iter ", dav_iter
        do run = 1, 6
            call run_smogd()
        end do
        write (*, "(a70,a,i4,a,i4)") trim(label), ":", n_fail_cfg, " failed of", n_tot_cfg
    end do
!
    write (*, "(/,a,i0,a,i0,/)") "non-symmetric and linear response: ", n_fail, " failures out of ", n_tot
contains
!
    subroutine build_nosym(deg_in, cplx_in)
        real(dp), intent(in) :: deg_in, cplx_in
        real(dp) :: u(n), w(n), d(n, n), fac
        call random_number(u); u = (u - 0.5_dp)*0.3_dp/sqrt(real(n, dp))
        call random_number(w); w = (w - 0.5_dp)*0.5_dp/sqrt(real(n, dp))
        d = 0.0_dp
        do i = 1, n
            d(i, i) = 1.0_dp + real(i - 1, dp)
            if (deg_in .gt. 0.0_dp .and. mod(i, 2_ip) .eq. 0) d(i, i) = d(i - 1, i - 1) + deg_in
        end do
        if (cplx_in .gt. 0.0_dp) then          ! a complex pair among the lowest eigenvalues
            d(1, 2) = cplx_in
            d(2, 1) = -cplx_in
        end if
        fac = 1.0_dp + dot_product(w, u)
        do j = 1, n
            do i = 1, n
                m_r(i, j) = d(i, j) + u(i)*dot_product(w, d(:, j)) &
                            - (d(i, 1)*0.0_dp)
            end do
        end do
!       m_r = (I + u w^T) D (I - u w^T / fac)
        do j = 1, n
            m_r(:, j) = m_r(:, j) - matmul(m_r, u)*w(j)/fac
        end do
        m_l = transpose(m_r)
    end subroutine build_nosym
!
    subroutine reference_nosym()
        real(dp) :: c(n, n)
        real(dp) :: wr_l(n), wi_l(n)
        c = m_r
        call dgeev('n', 'n', n, c, n, wr_l, wi_l, vdum, 1_ip, vdum, 1_ip, work, lwork, info)
        if (info .ne. 0) then
            write (*, *) "dgeev failed", info
            stop 1
        end if
!       keep the real eigenvalues, sorted
        n_kept = 0
        do i = 1, n
            if (abs(wi_l(i)) .lt. 1.0e-10_dp) then
                n_kept = n_kept + 1
                ref(n_kept) = wr_l(i)
            end if
        end do
        call sort(ref(1:n_kept))
    end subroutine reference_nosym
!
    subroutine run_nosym(side)
        character(len=2), intent(in) :: side
        if (allocated(eig)) deallocate (eig, evec, evec2)
        allocate (eig(n_max), evec(n, n_max), evec2(n, n_max))
        eig = 0.0_dp; evec = 0.0_dp; evec2 = 0.0_dp
        if (2*n_max .ge. n) return
        if (side .eq. "LR") then
            call dgl_davidson_nosym_driver(n, n_targ, n_max, mv_r, mv_l, pc_nosym, side, eig, evec, ok, &
                                           evec_2=evec2, dgl_tol=tol, dgl_max_iter=400_ip, &
                                           dgl_dav_iter=dav_iter, dgl_info=info)
        else
            call dgl_davidson_nosym_driver(n, n_targ, n_max, mv_r, mv_l, pc_nosym, side, eig, evec, ok, &
                                           dgl_tol=tol, dgl_max_iter=400_ip, dgl_dav_iter=dav_iter, dgl_info=info)
        end if
        err = maxval(abs(eig(1:n_targ) - ref(1:n_targ)))
        res = res_nosym(side)
        n_tot = n_tot + 1
        n_tot_cfg = n_tot_cfg + 1
        if (.not. ok .or. info .ne. dgl_success .or. err .gt. 1.0e-7_dp .or. res .gt. 1.0e-5_dp) then
            n_fail = n_fail + 1
            n_fail_cfg = n_fail_cfg + 1
            write (*, "(a,a,a,l2,a,i3,a,es9.2,a,es9.2)") "   FAIL nosym ", side, " ok", ok, " info", info, &
                " eig err", err, " residual", res
        end if
    end subroutine run_nosym
!
    real(dp) function res_nosym(side)
        character(len=2), intent(in) :: side
        real(dp) :: ax(n)
        integer(ip) :: kk
        res_nosym = 0.0_dp
        do kk = 1, n_targ
            if (side .eq. "L ") then
                ax = matmul(m_l, evec(:, kk))
            else
                ax = matmul(m_r, evec(:, kk))
            end if
            res_nosym = max(res_nosym, norm2(ax - eig(kk)*evec(:, kk))/norm2(evec(:, kk)))
            if (side .eq. "LR") then
                ax = matmul(m_l, evec2(:, kk))
                res_nosym = max(res_nosym, norm2(ax - eig(kk)*evec2(:, kk))/norm2(evec2(:, kk)))
            end if
        end do
    end function res_nosym
!
    subroutine build_lr(dpart)
        real(dp), intent(in) :: dpart
        real(dp) :: r
        integer(ip) :: ii, jj
        do jj = 1, n
            do ii = 1, n
                r = 1.0_dp/real(ii + jj, dp)
                apb(ii, jj) = r
                amb(ii, jj) = 0.2_dp*r
                spd(ii, jj) = merge(1.0_dp, 0.0_dp, ii .eq. jj)
                smd(ii, jj) = merge(1.0_dp, 0.0_dp, ii .eq. jj)
                if (dpart .gt. 0.0_dp .and. ii .ne. jj) then
                    spd(ii, jj) = spd(ii, jj) + sign(dpart/real(ii + jj, dp), real(jj - ii, dp))
                    smd(ii, jj) = smd(ii, jj) - sign(dpart/real(ii + jj, dp), real(jj - ii, dp))
                end if
            end do
            apb(jj, jj) = apb(jj, jj) + 5.0_dp + real(jj, dp)
            amb(jj, jj) = amb(jj, jj) + 2.0_dp + real(jj, dp)
        end do
    end subroutine build_lr
!
    subroutine reference_lr()
!       E z = omega M z with E = [A B; B A], M = [S D; -D -S]
        real(dp) :: aa(n, n), bb(n, n), ss(n, n), dd(n, n)
        integer(ip) :: ii, jj
        aa = 0.5_dp*(apb + amb)
        bb = 0.5_dp*(apb - amb)
        ss = 0.5_dp*(spd + smd)
        dd = 0.5_dp*(spd - smd)
        e2(1:n, 1:n) = aa; e2(1:n, n + 1:2*n) = bb
        e2(n + 1:2*n, 1:n) = bb; e2(n + 1:2*n, n + 1:2*n) = aa
        met(1:n, 1:n) = ss; met(1:n, n + 1:2*n) = dd
        met(n + 1:2*n, 1:n) = -dd; met(n + 1:2*n, n + 1:2*n) = -ss
        acopy = e2
        bcopy = met
        call dggev('n', 'n', 2*n, acopy, 2*n, bcopy, 2*n, wr, wi, ref_lr, vdum, 1_ip, vdum, 1_ip, work, lwork, info)
        if (info .ne. 0) then
            write (*, *) "dggev failed", info
            stop 1
        end if
!       eigenvalues are wr/ref_lr (beta); keep the positive real ones
        n_kept = 0
        do ii = 1, 2*n
            if (abs(wi(ii)) .lt. 1.0e-8_dp .and. abs(ref_lr(ii)) .gt. 1.0e-10_dp) then
                if (wr(ii)/ref_lr(ii) .gt. 1.0e-8_dp) then
                    n_kept = n_kept + 1
                    ref(n_kept) = wr(ii)/ref_lr(ii)
                end if
            end if
        end do
        call sort(ref(1:n_kept))
    end subroutine reference_lr
!
    subroutine run_smogd()
        if (allocated(evlr)) deallocate (evlr)
        if (allocated(eig)) deallocate (eig)
        allocate (eig(n_max), evlr(2*n, n_max))
        eig = 0.0_dp; evlr = 0.0_dp
        call dgl_smogd_driver(2*n, n_targ, n_max, mv_apb, mv_amb, mv_spd, mv_smd, pc_lr, eig, evlr, ok, &
                              dgl_tol=tol, dgl_max_iter=400_ip, dgl_dav_iter=dav_iter, dgl_info=info)
        err = maxval(abs(eig(1:n_targ) - ref(1:n_targ)))
        n_tot = n_tot + 1
        n_tot_cfg = n_tot_cfg + 1
        if (.not. ok .or. info .ne. dgl_success .or. err .gt. 1.0e-7_dp) then
            n_fail = n_fail + 1
            n_fail_cfg = n_fail_cfg + 1
            write (*, "(a,l2,a,i3,a,es9.2)") "   FAIL smogd ok", ok, " info", info, " eig err", err
        end if
    end subroutine run_smogd
!
    subroutine sort(x)
        real(dp), intent(inout) :: x(:)
        integer(ip) :: ii, jj
        real(dp) :: t
        do ii = 2, size(x, kind=ip)
            t = x(ii)
            jj = ii - 1
            do while (jj .ge. 1)
                if (x(jj) .le. t) exit
                x(jj + 1) = x(jj)
                jj = jj - 1
            end do
            x(jj + 1) = t
        end do
    end subroutine sort
!
    subroutine mv_r(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, m_r, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_r
    subroutine mv_l(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, m_l, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_l
    subroutine mv_apb(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, apb, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_apb
    subroutine mv_amb(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, amb, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_amb
    subroutine mv_spd(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, spd, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_spd
    subroutine mv_smd(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call dgemm('n', 'n', nn, m, nn, 1.0_dp, smd, n, x, nn, 0.0_dp, y, nn)
    end subroutine mv_smd
    subroutine pc_nosym(nn, m, shift, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: kk, ll
        real(dp) :: den
        do ll = 1, m
            do kk = 1, nn
                den = m_r(kk, kk) + shift
                if (abs(den) .lt. 1.0e-3_dp) den = sign(1.0e-3_dp, den)
                y(kk, ll) = x(kk, ll)/den
            end do
        end do
    end subroutine pc_nosym
    subroutine pc_lr(nn, m, fac, xp, xm, yp, ym)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: fac
        real(dp), intent(in) :: xp(nn, m), xm(nn, m)
        real(dp), intent(inout) :: yp(nn, m), ym(nn, m)
        integer(ip) :: kk, ll
        real(dp) :: den
        do ll = 1, m
            do kk = 1, nn
                den = fac*fac*apb(kk, kk)*amb(kk, kk) - 1.0_dp
                if (abs(den) .lt. 1.0e-3_dp) den = 1.0e-3_dp
                yp(kk, ll) = (fac*amb(kk, kk)*xp(kk, ll) + xm(kk, ll))/den
                ym(kk, ll) = (fac*apb(kk, kk)*xm(kk, ll) + xp(kk, ll))/den
            end do
        end do
    end subroutine pc_lr
!
end program stress_nosym
