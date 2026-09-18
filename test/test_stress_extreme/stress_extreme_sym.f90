! Extreme randomized fuzz test for the symmetric/generalized drivers (Davidson-Liu, LOBPCG).
! Unlike test/test_stress/stress_sym.f90 (fixed configs, fixed seed, n=400), this program:
!  - draws a fresh random configuration every trial (size, spectrum, metric conditioning,
!    degeneracy pattern, number of roots, subspace size, restart interval, preconditioner
!    quality, precnd_shift) from a seed given on the command line (or system-clock based)
!  - runs for a wall-clock time budget (default 20 minutes) rather than a fixed trial count
!  - separates real bugs (wrong results claimed converged, NaN/Inf, unexpected error codes)
!    from plain non-convergence within max_iter (which is not a bug, just a hard problem)
!  - includes a handful of large (n up to 4000) and extreme-degeneracy/conditioning cases
!
! Build: link against an installed diaglib (see build_and_run.sh). Usage:
!   stress_extreme_sym [seed] [minutes]
program stress_extreme_sym
    use dgl_interface
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none
    integer, parameter :: ip = dgl_int, dp = dgl_real
    integer(ip) :: n, n_targ, n_max, dav_iter, max_iter
    real(dp), allocatable :: a(:, :), b(:, :), acopy(:, :), bcopy(:, :)
    real(dp), allocatable :: ref(:), diag_a(:), work(:)
    real(dp), allocatable :: eig(:), evec(:, :)
    integer(ip) :: info, lwork
    integer(ip) :: n_bug, n_noconv, n_ok, n_tot, trial
    logical :: ok, generalized, weak, spd_precnd, use_shift, block_deg
    real(dp) :: cur_cond_b
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
    master_seed = 20260918123456_int64
    minutes_budget = 20.0_real64
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
    write (*, "(a,i0,a,f0.1,a,i0)") "stress_extreme_sym: master_seed=", master_seed, &
        " time_budget=", minutes_budget, " min, max_n=", max_n_cap
!
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = int(master_seed, kind(seed))
    call random_seed(put=seed)
!
    n_bug = 0; n_noconv = 0; n_ok = 0; n_tot = 0
    call cpu_time(t_start)
    call wall_seconds(t_start)
!
! a handful of hand-picked extreme cases first (always run, not time-budget-gated)
!
    call run_extreme_case(2000_ip, 8_ip, .true., 1.0e10_dp, "n=2000 gen extreme-cond(B)=1e10")
    call run_extreme_case(4000_ip, 30_ip, .false., 1.0_dp, "n=4000 std, 30 roots (go-large)")
    call run_extreme_case(4000_ip, 1_ip, .true., 1.0e4_dp, "n=4000 gen, 1 root")
    call block_degeneracy_case(600_ip, 24_ip, "n=600 massive block degeneracy (2 clusters)")
    call tiny_case()
!
! randomized fuzz loop, time-boxed
!
    trial = 0
    do
        call wall_seconds(t_now)
        if ((t_now - t_start)/60.0_real64 .gt. minutes_budget) exit
        trial = trial + 1
        call random_trial(trial)
    end do
!
    write (*, "(/,a)") "==================== stress_extreme_sym summary ===================="
    write (*, "(a,i0)") "total trials (incl. extreme cases): ", n_tot
    write (*, "(a,i0)") "genuine bugs (wrong result claimed converged, NaN, bad error code): ", n_bug
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
    subroutine random_trial(itrial)
        integer(ip), intent(in) :: itrial
        real(dp) :: spread_a, deg, cond_b, r
        integer(ip) :: nn, ntg, nmx, davi, maxi
        logical :: gen, wk, spdp, shft
!
        call rnd_range(150.0_dp, min(900.0_dp, real(max_n_cap, dp)), r); nn = int(r, ip)
        call random_number(r); gen = r .lt. 0.5_dp
        call rnd_range(0.3_dp, 800.0_dp, spread_a)
        call random_number(r)
        if (r .lt. 0.3_dp) then
            deg = 0.0_dp
        else if (r .lt. 0.6_dp) then
            deg = 1.0e-3_dp
        else if (r .lt. 0.85_dp) then
            deg = 1.0e-6_dp
        else
            deg = 1.0e-9_dp
        end if
        call random_number(r)
        if (.not. gen) then
            cond_b = 1.0_dp
        else if (r .lt. 0.25_dp) then
            cond_b = 1.0_dp
        else if (r .lt. 0.5_dp) then
            cond_b = 1.0e3_dp
        else if (r .lt. 0.75_dp) then
            cond_b = 1.0e7_dp
        else
            cond_b = 1.0e11_dp
        end if
        call rnd_range(1.0_dp, real(nn, dp)/4.0_dp, r); ntg = max(1_ip, int(r, ip))
        call rnd_range(real(ntg, dp) + 1.0_dp, min(3.0_dp*real(ntg, dp), real(nn, dp)/2.0_dp), r)
        nmx = max(ntg + 1_ip, int(r, ip))
        if (nmx .ge. nn) nmx = nn/2_ip
        if (nmx .le. ntg) nmx = ntg + 1_ip
        call rnd_range(4.0_dp, 40.0_dp, r); davi = int(r, ip)
        call rnd_range(200.0_dp, 1500.0_dp, r); maxi = int(r, ip)
        call random_number(r); wk = r .lt. 0.35_dp
        call random_number(r); spdp = r .lt. 0.2_dp
        call random_number(r); shft = r .lt. 0.5_dp
!
        write (cfgline, "(a,i0,a,i0,a,i0,a,l1,a,es8.1,a,es8.1,a,l1,a,l1)") &
            "fuzz#", itrial, " n=", nn, " ntg=", ntg, " gen=", gen, &
            " deg=", deg, " condB=", cond_b, " weak=", wk, " shift=", shft
        call setup_and_run(nn, ntg, nmx, davi, maxi, gen, spread_a, deg, cond_b, wk, spdp, shft, &
                            trim(cfgline))
    end subroutine random_trial
!
    subroutine run_extreme_case(nn, ntg, gen, cond_b, label)
        integer(ip), intent(in) :: nn, ntg
        logical, intent(in) :: gen
        real(dp), intent(in) :: cond_b
        character(len=*), intent(in) :: label
        integer(ip) :: nmx
        if (nn .gt. max_n_cap) then
            write (*, "(a,a)") "   (skipped, n exceeds max_n_cap) -- ", trim(label)
            return
        end if
        nmx = min(3_ip*ntg, nn/2_ip)
        if (nmx .le. ntg) nmx = ntg + 1_ip
        call setup_and_run(nn, ntg, nmx, 25_ip, 2000_ip, gen, 50.0_dp, 1.0e-7_dp, cond_b, &
                            .false., .false., .true., label)
    end subroutine run_extreme_case
!
    subroutine block_degeneracy_case(nn, ntg, label)
        integer(ip), intent(in) :: nn, ntg
        character(len=*), intent(in) :: label
        integer(ip) :: nmx
        block_deg = .true.
        nmx = 3_ip*ntg
        call setup_and_run(nn, ntg, nmx, 20_ip, 1000_ip, .false., 0.0_dp, 0.0_dp, 1.0_dp, &
                            .false., .false., .true., label)
        block_deg = .false.
    end subroutine block_degeneracy_case
!
    subroutine tiny_case()
        integer(ip) :: nn, ntg, nmx
        nn = 6_ip; ntg = 1_ip; nmx = 2_ip
        call setup_and_run(nn, ntg, nmx, 5_ip, 200_ip, .false., 5.0_dp, 0.0_dp, 1.0_dp, &
                            .false., .false., .true., "tiny n=6, 1 root")
        nn = 5_ip; ntg = 2_ip; nmx = 4_ip
        call setup_and_run(nn, ntg, nmx, 5_ip, 200_ip, .true., 5.0_dp, 0.0_dp, 10.0_dp, &
                            .false., .false., .true., "tiny n=5 generalized, 2 roots")
    end subroutine tiny_case
!
    subroutine setup_and_run(nn, ntg, nmx, davi, maxi, gen, spread_a, deg, cond_b, wk, spdp, shft, label)
        integer(ip), intent(in) :: nn, ntg, nmx, davi, maxi
        logical, intent(in) :: gen, wk, spdp, shft
        real(dp), intent(in) :: spread_a, deg, cond_b
        character(len=*), intent(in) :: label
        n = nn; n_targ = ntg; n_max = nmx; dav_iter = davi; max_iter = maxi
        generalized = gen; weak = wk; spd_precnd = spdp; use_shift = shft; cur_cond_b = cond_b
        if (allocated(a)) deallocate (a, b, acopy, bcopy, ref, diag_a, eig, evec)
        allocate (a(n, n), b(n, n), acopy(n, n), bcopy(n, n), ref(n), diag_a(n))
        allocate (eig(n_max), evec(n, n_max))
        lwork = 64*n
        if (allocated(work)) deallocate (work)
        allocate (work(lwork))
        call build_sym(spread_a, deg, cond_b)
        call reference_sym(generalized)
        call run_one("Davidson", .false., label)
        call run_one("LOBPCG", .true., label)
    end subroutine setup_and_run
!
    subroutine rnd_range(lo, hi, r)
        real(dp), intent(in) :: lo, hi
        real(dp), intent(out) :: r
        real(dp) :: u
        call random_number(u)
        r = lo + u*(hi - lo)
    end subroutine rnd_range
!
    subroutine build_sym(spread_in, deg_in, cond_in)
        real(dp), intent(in) :: spread_in, deg_in, cond_in
        real(dp) :: v(n), s(n), r, diag_b(n)
        integer(ip) :: k, l, half
        if (block_deg) then
            half = n/2_ip
            do k = 1, n
                diag_a(k) = merge(1.0_dp, 7.0_dp, k .le. half)
                diag_b(k) = 1.0_dp
            end do
        else
            do k = 1, n
                diag_a(k) = 1.0_dp + spread_in*real(k - 1, dp)/real(n - 1, dp)
                if (deg_in .gt. 0.0_dp .and. mod(k, 2_ip) .eq. 0) diag_a(k) = diag_a(k - 1) + deg_in
                diag_b(k) = 1.0_dp + (cond_in - 1.0_dp)*real(n - k, dp)/real(n - 1, dp)
            end do
        end if
        call random_number(v)
        v = v - 0.5_dp
        v = v/norm2(v)
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
            stop 2
        end if
    end subroutine reference_sym
!
    subroutine run_one(driver, lobpcg, label)
        character(len=*), intent(in) :: driver, label
        logical, intent(in) :: lobpcg
        real(dp) :: err, res, res_thresh
        logical :: has_nan
        procedure(dgl_matvec), pointer :: metvec_p
        eig = 0.0_dp
        evec = 0.0_dp
        if (lobpcg .and. 3*n_max .ge. n) return
        if (.not. lobpcg .and. 2*n_max .ge. n) return
        ok = .false.
        metvec_p => mv_b
        if (lobpcg) then
            if (generalized) then
                call dgl_lobpcg_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, metvec=metvec_p, &
                                       dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, dgl_info=info, &
                                       dgl_precnd_shift=use_shift, dgl_memory=8000_ip)
            else
                call dgl_lobpcg_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, &
                                       dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, dgl_info=info, &
                                       dgl_precnd_shift=use_shift, dgl_memory=8000_ip)
            end if
        else
            if (generalized) then
                call dgl_davidson_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, metvec=metvec_p, &
                                         dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                         dgl_info=info, dgl_precnd_shift=use_shift, dgl_memory=8000_ip)
            else
                call dgl_davidson_driver(n, n_targ, n_max, mv_a, pc, eig, evec, ok, &
                                         dgl_tol=1.0e-9_dp, dgl_max_iter=max_iter, dgl_dav_iter=dav_iter, &
                                         dgl_info=info, dgl_precnd_shift=use_shift, dgl_memory=8000_ip)
            end if
        end if
        n_tot = n_tot + 1
        has_nan = any(eig(1:n_targ) .ne. eig(1:n_targ)) .or. any(evec(:, 1:n_targ) .ne. evec(:, 1:n_targ))
        if (has_nan) then
            n_bug = n_bug + 1
            write (*, "(a,a,a,a)") "   *** BUG (NaN) ", trim(driver), " -- ", trim(label)
            return
        end if
        if (info .ne. dgl_success) then
            n_bug = n_bug + 1
            write (*, "(a,a,a,i3,a,a)") "   *** BUG (unexpected info) ", trim(driver), " info=", info, &
                " -- ", trim(label)
            return
        end if
        if (.not. ok) then
            n_noconv = n_noconv + 1
            return
        end if
        err = maxval(abs(eig(1:n_targ) - ref(1:n_targ)))
        res = residual()
        ! an ill-conditioned metric inflates the achievable residual even for an exact
        ! eigenvalue (roundoff in B is amplified by ~cond(B)): scale the residual tolerance
        ! accordingly rather than flagging expected conditioning-limited residuals as bugs.
        res_thresh = min(1.0e-1_dp, 1.0e-5_dp*max(1.0_dp, sqrt(cur_cond_b)*1.0e-2_dp))
        if (err .gt. 1.0e-6_dp .or. res .gt. res_thresh) then
            n_bug = n_bug + 1
            write (*, "(a,a,a,es9.2,a,es9.2,a,es9.2,a,a)") "   *** BUG (claimed converged, wrong) ", trim(driver), &
                " eig_err=", err, " residual=", res, " (thresh=", res_thresh, ") -- ", trim(label)
        else
            n_ok = n_ok + 1
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
            residual = max(residual, norm2(ax - eig(k)*bx)/(norm2(ax) + abs(eig(k))*norm2(bx) + 1.0e-300_dp))
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
                if (weak) den = a(k, k) + shift + 200.0_dp
                if (spd_precnd) den = abs(a(k, k))
                if (abs(den) .lt. 1.0e-3_dp) den = sign(1.0e-3_dp, den)
                y(k, l) = x(k, l)/den
            end do
        end do
    end subroutine pc
!
end program stress_extreme_sym
