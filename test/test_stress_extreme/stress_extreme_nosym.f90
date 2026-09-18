! Extreme randomized fuzz test for the non-symmetric Davidson and SMO-GD drivers.
! Same philosophy as stress_extreme_sym.f90: randomized configs each trial, time-boxed,
! plus a few hand-picked extreme cases (larger n, heavier non-normality, many complex pairs,
! extreme conditioning), and a clear separation of genuine bugs from plain non-convergence.
program stress_extreme_nosym
    use dgl_interface
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none
    integer, parameter :: ip = dgl_int, dp = dgl_real
    integer(ip) :: n
    real(dp), allocatable :: m_r(:, :), m_l(:, :), apb(:, :), amb(:, :), spd(:, :), smd(:, :)
    real(dp), allocatable :: e2(:, :), met(:, :), acopy(:, :), bcopy(:, :)
    real(dp), allocatable :: ref(:), ref_lr(:), wr(:), wi(:), work(:), vdum(:, :)
    real(dp), allocatable :: eig(:), evec(:, :), evec2(:, :)
    integer(ip) :: info, lwork, n_targ, n_max, dav_iter, max_iter, n_kept
    integer(ip) :: n_bug, n_noconv, n_ok, n_tot, trial
    logical :: ok
    character(len=200) :: cfgline
    integer :: seed_size
    integer, allocatable :: seed(:)
    integer(int64) :: master_seed
    real(real64) :: t_start, t_now, minutes_budget
    integer :: nargs, ios
    character(len=64) :: argbuf
    integer(ip) :: max_n_cap
!
    nargs = command_argument_count()
    master_seed = 987654321_int64
    minutes_budget = 15.0_real64
    max_n_cap = 100000_ip
    if (nargs .ge. 1) then
        call get_command_argument(1, argbuf)
        read (argbuf, *, iostat=ios) master_seed
    end if
    if (nargs .ge. 2) then
        call get_command_argument(2, argbuf)
        read (argbuf, *, iostat=ios) minutes_budget
    end if
    if (nargs .ge. 3) then
        call get_command_argument(3, argbuf)
        read (argbuf, *, iostat=ios) max_n_cap
    end if
    write (*, "(a,i0,a,f0.1,a,i0)") "stress_extreme_nosym: master_seed=", master_seed, &
        " time_budget=", minutes_budget, " min, max_n=", max_n_cap
!
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = int(master_seed, kind(seed))
    call random_seed(put=seed)
!
    n_bug = 0; n_noconv = 0; n_ok = 0; n_tot = 0
    call wall_seconds(t_start)
!
    call extreme_nosym_case(1500_ip, 20_ip, 1.0e-6_dp, 0.6_dp, "n=1500 nosym, heavy complex content")
    call extreme_nosym_case(1500_ip, 1_ip, 0.0_dp, 0.0_dp, "n=1500 nosym, 1 root")
    call extreme_lr_case(800_ip, 15_ip, 0.35_dp, "n2=1600 smogd, strong D antisymmetric part")
!
    trial = 0
    do
        call wall_seconds(t_now)
        if ((t_now - t_start)/60.0_real64 .gt. minutes_budget) exit
        trial = trial + 1
        call random_nosym_trial(trial)
        call wall_seconds(t_now)
        if ((t_now - t_start)/60.0_real64 .gt. minutes_budget) exit
        call random_lr_trial(trial)
    end do
!
    write (*, "(/,a)") "================== stress_extreme_nosym summary ===================="
    write (*, "(a,i0)") "total trials (incl. extreme cases): ", n_tot
    write (*, "(a,i0)") "genuine bugs: ", n_bug
    write (*, "(a,i0)") "plain non-convergence (not a bug): ", n_noconv
    write (*, "(a,i0)") "converged and correct: ", n_ok
    write (*, "(a,i0)") "random fuzz trials run: ", trial
    if (n_bug .gt. 0) then
        write (*, "(a)") "RESULT: BUGS FOUND"
        stop 1
    else
        write (*, "(a)") "RESULT: CLEAN"
    end if
!
contains
!
    subroutine wall_seconds(t)
        real(real64), intent(out) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, real64)/real(rate, real64)
    end subroutine wall_seconds
!
    subroutine rnd_range(lo, hi, r)
        real(dp), intent(in) :: lo, hi
        real(dp), intent(out) :: r
        real(dp) :: u
        call random_number(u)
        r = lo + u*(hi - lo)
    end subroutine rnd_range
!
    subroutine random_nosym_trial(itrial)
        integer(ip), intent(in) :: itrial
        real(dp) :: deg, cplx, r
        integer(ip) :: nn, ntg, nmx, davi, maxi, iside
        character(len=2) :: sides(3) = ["R ", "L ", "LR"]
        call rnd_range(80.0_dp, min(500.0_dp, real(max_n_cap, dp)), r); nn = int(r, ip)
        call random_number(r)
        deg = merge(0.0_dp, 1.0e-5_dp, r .lt. 0.4_dp)
        call random_number(r)
        cplx = merge(0.0_dp, merge(0.05_dp, 0.7_dp, r .lt. 0.7_dp), r .lt. 0.4_dp)
        call rnd_range(1.0_dp, real(nn, dp)/6.0_dp, r); ntg = max(1_ip, int(r, ip))
        call rnd_range(real(ntg, dp) + 1.0_dp, min(3.0_dp*real(ntg, dp), real(nn, dp)/3.0_dp), r)
        nmx = max(ntg + 1_ip, int(r, ip))
        if (2*nmx .ge. nn) nmx = nn/3_ip
        if (nmx .le. ntg) nmx = ntg + 1_ip
        call rnd_range(4.0_dp, 35.0_dp, r); davi = int(r, ip)
        call rnd_range(200.0_dp, 1200.0_dp, r); maxi = int(r, ip)
        call random_number(r); iside = 1 + int(r*3.0_dp)
        if (iside .gt. 3) iside = 3
        write (cfgline, "(a,i0,a,i0,a,i0,a,es8.1,a,es8.1,a,a)") &
            "nosym-fuzz#", itrial, " n=", nn, " ntg=", ntg, " deg=", deg, " cplx=", cplx, &
            " side=", trim(sides(iside))
        call setup_nosym(nn, deg, cplx)
        call run_nosym_case(ntg, nmx, davi, maxi, sides(iside), trim(cfgline))
    end subroutine random_nosym_trial
!
    subroutine random_lr_trial(itrial)
        integer(ip), intent(in) :: itrial
        real(dp) :: dpart, r
        integer(ip) :: nn, ntg, nmx, davi, maxi
        call rnd_range(60.0_dp, min(400.0_dp, real(max_n_cap, dp)), r); nn = int(r, ip)
        call rnd_range(0.0_dp, 0.45_dp, dpart)
        call rnd_range(1.0_dp, real(nn, dp)/6.0_dp, r); ntg = max(1_ip, int(r, ip))
        call rnd_range(real(ntg, dp) + 1.0_dp, min(3.0_dp*real(ntg, dp), real(nn, dp)/3.0_dp), r)
        nmx = max(ntg + 1_ip, int(r, ip))
        call rnd_range(4.0_dp, 35.0_dp, r); davi = int(r, ip)
        call rnd_range(200.0_dp, 1200.0_dp, r); maxi = int(r, ip)
        write (cfgline, "(a,i0,a,i0,a,i0,a,es8.1)") &
            "smogd-fuzz#", itrial, " n2=", 2*nn, " ntg=", ntg, " dpart=", dpart
        call setup_lr(nn, dpart)
        call run_lr_case(ntg, nmx, davi, maxi, trim(cfgline))
    end subroutine random_lr_trial
!
    subroutine extreme_nosym_case(nn, ntg, deg, cplx, label)
        integer(ip), intent(in) :: nn, ntg
        real(dp), intent(in) :: deg, cplx
        character(len=*), intent(in) :: label
        integer(ip) :: nmx
        if (nn .gt. max_n_cap) then
            write (*, "(a,a)") "   (skipped, n exceeds max_n_cap) -- ", trim(label)
            return
        end if
        nmx = min(3_ip*ntg, nn/3_ip)
        if (nmx .le. ntg) nmx = ntg + 1_ip
        call setup_nosym(nn, deg, cplx)
        call run_nosym_case(ntg, nmx, 25_ip, 2000_ip, "LR", label)
    end subroutine extreme_nosym_case
!
    subroutine extreme_lr_case(nn, ntg, dpart, label)
        integer(ip), intent(in) :: nn, ntg
        real(dp), intent(in) :: dpart
        character(len=*), intent(in) :: label
        integer(ip) :: nmx
        if (2_ip*nn .gt. max_n_cap) then
            write (*, "(a,a)") "   (skipped, n exceeds max_n_cap) -- ", trim(label)
            return
        end if
        nmx = 3_ip*ntg
        call setup_lr(nn, dpart)
        call run_lr_case(ntg, nmx, 25_ip, 2000_ip, label)
    end subroutine extreme_lr_case
!
    subroutine setup_nosym(nn, deg_in, cplx_in)
        integer(ip), intent(in) :: nn
        real(dp), intent(in) :: deg_in, cplx_in
        real(dp) :: fac
        real(dp), allocatable :: u(:), w(:), d(:, :)
        integer(ip) :: i, j
        n = nn
        if (allocated(m_r)) deallocate (m_r, m_l)
        allocate (m_r(n, n), m_l(n, n))
        allocate (u(n), w(n), d(n, n))
        call random_number(u); u = (u - 0.5_dp)*0.3_dp/sqrt(real(n, dp))
        call random_number(w); w = (w - 0.5_dp)*0.5_dp/sqrt(real(n, dp))
        d = 0.0_dp
        do i = 1, n
            d(i, i) = 1.0_dp + real(i - 1, dp)
            if (deg_in .gt. 0.0_dp .and. mod(i, 2) .eq. 0) d(i, i) = d(i - 1, i - 1) + deg_in
        end do
        if (cplx_in .gt. 0.0_dp) then
            ! complex pairs only among the lowest few eigenvalues (like the repo's own
            ! stress_nosym.f90): injecting them across the whole spectrum with unit diagonal
            ! spacing hits an exact discriminant-zero coincidence at cplx_in=0.5 (turning the
            ! "complex" pairs into massively degenerate real ones instead) and made the dense
            ! LAPACK reference's real/complex split disagree with which roots the driver tracks
            do i = 1, min(6_ip, n - 1_ip), 2
                d(i, i + 1) = cplx_in
                d(i + 1, i) = -cplx_in
            end do
        end if
        fac = 1.0_dp + dot_product(w, u)
        do j = 1, n
            do i = 1, n
                m_r(i, j) = d(i, j) + u(i)*dot_product(w, d(:, j))
            end do
        end do
        do j = 1, n
            m_r(:, j) = m_r(:, j) - matmul(m_r, u)*w(j)/fac
        end do
        m_l = transpose(m_r)
        call setup_workspace(n)
        call reference_nosym()
    end subroutine setup_nosym
!
    subroutine setup_workspace(nn)
        integer(ip), intent(in) :: nn
        lwork = 64*2*nn
        if (allocated(work)) deallocate (work)
        allocate (work(lwork))
        if (allocated(vdum)) deallocate (vdum)
        allocate (vdum(1, 1))
        if (allocated(wr)) deallocate (wr, wi)
        allocate (wr(2*nn), wi(2*nn))
        if (allocated(ref)) deallocate (ref)
        allocate (ref(nn))
    end subroutine setup_workspace
!
    subroutine reference_nosym()
        real(dp), allocatable :: c(:, :), wr_l(:), wi_l(:)
        integer(ip) :: i
        allocate (c(n, n), wr_l(n), wi_l(n))
        c = m_r
        call dgeev('n', 'n', n, c, n, wr_l, wi_l, vdum, 1_ip, vdum, 1_ip, work, lwork, info)
        if (info .ne. 0) then
            write (*, *) "dgeev failed", info
            stop 2
        end if
        n_kept = 0
        do i = 1, n
            if (abs(wi_l(i)) .lt. 1.0e-10_dp) then
                n_kept = n_kept + 1
                ref(n_kept) = wr_l(i)
            end if
        end do
        call sort(ref(1:n_kept))
        deallocate (c, wr_l, wi_l)
    end subroutine reference_nosym
!
    subroutine run_nosym_case(ntg, nmx, davi, maxi, side, label)
        integer(ip), intent(in) :: ntg, nmx, davi, maxi
        character(len=2), intent(in) :: side
        character(len=*), intent(in) :: label
        real(dp) :: err, res
        logical :: has_nan
        n_targ = ntg; n_max = nmx; dav_iter = davi; max_iter = maxi
        if (ntg .gt. n_kept) return
        if (2*n_max .ge. n) return
        if (allocated(eig)) deallocate (eig, evec, evec2)
        allocate (eig(n_max), evec(n, n_max), evec2(n, n_max))
        eig = 0.0_dp; evec = 0.0_dp; evec2 = 0.0_dp
        ok = .false.
        if (side .eq. "LR") then
            call dgl_davidson_nosym_driver(n, n_targ, n_max, mv_r, mv_l, pc_nosym, side, eig, evec, ok, &
                                           evec_2=evec2, dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, &
                                           dgl_dav_iter=dav_iter, dgl_info=info, dgl_memory=8000_ip)
        else
            call dgl_davidson_nosym_driver(n, n_targ, n_max, mv_r, mv_l, pc_nosym, side, eig, evec, ok, &
                                           dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                           dgl_info=info, dgl_memory=8000_ip)
        end if
        n_tot = n_tot + 1
        has_nan = any(eig(1:n_targ) .ne. eig(1:n_targ))
        if (has_nan) then
            n_bug = n_bug + 1
            write (*, "(a,a,a,a)") "   *** BUG (NaN) nosym ", side, " -- ", trim(label)
            return
        end if
        if (info .ne. dgl_success) then
            n_bug = n_bug + 1
            write (*, "(a,i3,a,a,a,a)") "   *** BUG (unexpected info) nosym info=", info, " side=", &
                side, " -- ", trim(label)
            return
        end if
        if (.not. ok) then
            n_noconv = n_noconv + 1
            return
        end if
        err = maxval(abs(eig(1:n_targ) - ref(1:n_targ)))
        res = res_nosym(side)
        if (err .gt. 1.0e-6_dp .or. res .gt. 1.0e-5_dp) then
            n_bug = n_bug + 1
            write (*, "(a,a,a,es9.2,a,es9.2,a,a)") "   *** BUG (claimed converged, wrong) nosym ", side, &
                " eig_err=", err, " residual=", res, " -- ", trim(label)
        else
            n_ok = n_ok + 1
        end if
    end subroutine run_nosym_case
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
            res_nosym = max(res_nosym, norm2(ax - eig(kk)*evec(:, kk))/(norm2(evec(:, kk)) + 1.0e-300_dp))
            if (side .eq. "LR") then
                ax = matmul(m_l, evec2(:, kk))
                res_nosym = max(res_nosym, norm2(ax - eig(kk)*evec2(:, kk))/(norm2(evec2(:, kk)) + 1.0e-300_dp))
            end if
        end do
    end function res_nosym
!
    subroutine setup_lr(nn, dpart)
        integer(ip), intent(in) :: nn
        real(dp), intent(in) :: dpart
        real(dp) :: r
        integer(ip) :: ii, jj
        n = nn
        if (allocated(apb)) deallocate (apb, amb, spd, smd)
        allocate (apb(n, n), amb(n, n), spd(n, n), smd(n, n))
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
        if (allocated(e2)) deallocate (e2, met, acopy, bcopy)
        allocate (e2(2*n, 2*n), met(2*n, 2*n), acopy(2*n, 2*n), bcopy(2*n, 2*n))
        call setup_workspace(2*n)
        if (allocated(ref_lr)) deallocate (ref_lr)
        allocate (ref_lr(2*n))
        call reference_lr()
    end subroutine setup_lr
!
    subroutine reference_lr()
        real(dp), allocatable :: aa(:, :), bb(:, :), ss(:, :), dd(:, :)
        integer(ip) :: ii
        allocate (aa(n, n), bb(n, n), ss(n, n), dd(n, n))
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
            stop 2
        end if
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
        deallocate (aa, bb, ss, dd)
    end subroutine reference_lr
!
    subroutine run_lr_case(ntg, nmx, davi, maxi, label)
        integer(ip), intent(in) :: ntg, nmx, davi, maxi
        character(len=*), intent(in) :: label
        real(dp) :: err
        logical :: has_nan
        n_targ = ntg; n_max = nmx; dav_iter = davi; max_iter = maxi
        if (ntg .gt. n_kept) return
        if (allocated(eig)) deallocate (eig)
        if (allocated(evec)) deallocate (evec)
        allocate (eig(n_max), evec(2*n, n_max))
        eig = 0.0_dp; evec = 0.0_dp
        ok = .false.
        call dgl_smogd_driver(2*n, n_targ, n_max, mv_apb, mv_amb, mv_spd, mv_smd, pc_lr, eig, evec, ok, &
                              dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, dgl_info=info, &
                              dgl_memory=8000_ip)
        n_tot = n_tot + 1
        has_nan = any(eig(1:n_targ) .ne. eig(1:n_targ))
        if (has_nan) then
            n_bug = n_bug + 1
            write (*, "(a,a)") "   *** BUG (NaN) smogd -- ", trim(label)
            return
        end if
        if (info .ne. dgl_success) then
            n_bug = n_bug + 1
            write (*, "(a,i3,a,a)") "   *** BUG (unexpected info) smogd info=", info, " -- ", trim(label)
            return
        end if
        if (.not. ok) then
            n_noconv = n_noconv + 1
            return
        end if
        err = maxval(abs(eig(1:n_targ) - ref(1:n_targ)))
        if (err .gt. 1.0e-6_dp) then
            n_bug = n_bug + 1
            write (*, "(a,es9.2,a,a)") "   *** BUG (claimed converged, wrong) smogd eig_err=", err, &
                " -- ", trim(label)
        else
            n_ok = n_ok + 1
        end if
    end subroutine run_lr_case
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
end program stress_extreme_nosym
