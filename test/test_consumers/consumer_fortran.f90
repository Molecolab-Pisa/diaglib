! Fortran program using an installed DiagLib, as a user would: only "use dgl_interface", and
! the installed CMake package (find_package(diaglib)).
!
! Test problems with known eigenvalues, applied matrix-free:
!   symmetric      A = H diag(d) H,  H = I - 2 v v^T (Householder)       -> lambda = d
!   generalized    B = H diag(s) H                                       -> lambda = d/s
!   non-symmetric  M = P diag(d) P^-1, P = I + u w^T                      -> lambda = d
!   linear resp.   A+B = H diag(a) H, A-B = H diag(b) H, S = I, D = 0    -> omega = sqrt(a b)
module problems
    use dgl_interface, only: ip => dgl_int, dp => dgl_real
    implicit none
    integer(ip), parameter :: n = 1000
    real(dp) :: v(n), u(n), w(n), d(n), s(n), a(n), b(n), wu
    ! scale factor of the matrix, different for each thread in the concurrent test
    real(dp) :: scale = 1.0_dp
    integer :: nested_calls = 0
    !$omp threadprivate(scale)
contains
    subroutine setup()
        integer(ip) :: i
        do i = 1, n
            v(i) = sin(real(i, dp))
            u(i) = 0.3_dp*cos(real(2*i, dp))/sqrt(real(n, dp))
            w(i) = 0.5_dp*sin(real(3*i, dp))/sqrt(real(n, dp))
            d(i) = 1.0_dp + real(i - 1, dp) + 0.5_dp*sin(real(i, dp))**2
            s(i) = 1.0_dp + 0.2_dp*cos(real(i, dp))**2
            a(i) = 2.0_dp + real(i, dp)
            b(i) = 1.0_dp + 0.5_dp*real(i, dp)
        end do
        v = v/norm2(v)
        wu = dot_product(w, u)
    end subroutine setup
    subroutine house(m, x)          ! x = H x, column by column
        integer(ip), intent(in) :: m
        real(dp), intent(inout) :: x(n, m)
        integer(ip) :: j
        do j = 1, m
            x(:, j) = x(:, j) - 2.0_dp*dot_product(v, x(:, j))*v
        end do
    end subroutine house
    subroutine hdh(diag, nn, m, x, y)
        real(dp), intent(in) :: diag(n)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: j
        y = x
        call house(m, y)
        do j = 1, m
            y(:, j) = diag*y(:, j)
        end do
        call house(m, y)
    end subroutine hdh
    subroutine mv_sym(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call hdh(d, nn, m, x, y)
        y = scale*y
    end subroutine mv_sym
    subroutine mv_metric(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call hdh(s, nn, m, x, y)
    end subroutine mv_metric
    subroutine pc_sym(nn, m, shift, x, y)  ! approximate diagonal preconditioner
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: j
        real(dp) :: den(nn)
        den = scale*d + shift
        where (abs(den) .lt. 1.0e-3_dp) den = sign(1.0e-3_dp, den)
        do j = 1, m
            y(:, j) = x(:, j)/den
        end do
    end subroutine pc_sym
    subroutine pc_spd(nn, m, shift, x, y)   ! LOBPCG needs a positive definite preconditioner
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: shift
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: j
        if (.false.) y = shift
        do j = 1, m
            y(:, j) = x(:, j)/(scale*d)
        end do
    end subroutine pc_spd
    subroutine mv_r(nn, m, x, y)         ! P D P^-1 x, P^-1 = I - u w^T/(1 + w.u)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: j
        do j = 1, m
            y(:, j) = x(:, j) - u*dot_product(w, x(:, j))/(1.0_dp + wu)
            y(:, j) = d*y(:, j)
            y(:, j) = y(:, j) + u*dot_product(w, y(:, j))
        end do
    end subroutine mv_r
    subroutine mv_l(nn, m, x, y)         ! (P D P^-1)^T x = P^-T D P^T x
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        integer(ip) :: j
        do j = 1, m
            y(:, j) = x(:, j) + w*dot_product(u, x(:, j))
            y(:, j) = d*y(:, j)
            y(:, j) = y(:, j) - w*dot_product(u, y(:, j))/(1.0_dp + wu)
        end do
    end subroutine mv_l
    subroutine apb(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call hdh(a, nn, m, x, y)
    end subroutine apb
    subroutine amb(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        call hdh(b, nn, m, x, y)
    end subroutine amb
    subroutine ident(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        y = x
    end subroutine ident
    subroutine zero_mv(nn, m, x, y)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        if (.false.) y = x
        y = 0.0_dp
    end subroutine zero_mv
    subroutine lrprec(nn, m, fac, xp, xm, yp, ym)
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: fac
        real(dp), intent(in) :: xp(nn, m), xm(nn, m)
        real(dp), intent(inout) :: yp(nn, m), ym(nn, m)
        integer(ip) :: j
        real(dp) :: den(nn)
        den = fac*fac*a*b - 1.0_dp       ! crude diagonal model in the original basis
        where (abs(den) .lt. 1.0e-3_dp) den = 1.0e-3_dp
        do j = 1, m
            yp(:, j) = (fac*b*xp(:, j) + xm(:, j))/den
            ym(:, j) = (fac*a*xm(:, j) + xp(:, j))/den
        end do
    end subroutine lrprec
    subroutine mv_nested(nn, m, x, y)   ! calls DiagLib again before applying the matrix
        use dgl_interface
        integer(ip), intent(in) :: nn, m
        real(dp), intent(in) :: x(nn, m)
        real(dp), intent(inout) :: y(nn, m)
        real(dp) :: e(4), ev(n, 4)
        logical :: ok
        integer(ip) :: info
        ev = 0.0_dp
        !$omp atomic
        nested_calls = nested_calls + 1
        call dgl_lobpcg_driver(n, 2_ip, 4_ip, mv_sym, pc_spd, e, ev, ok, dgl_tol=1.0e-8_dp, dgl_info=info)
        if (.not. ok .or. info .ne. dgl_success .or. abs(e(1) - scale*minval(d)) .gt. 1.0e-8_dp) then
            error stop "nested call failed"
        end if
        call mv_sym(nn, m, x, y)
    end subroutine mv_nested
end module problems

program ext_fortran
    use dgl_interface
    use problems
    use omp_lib
    implicit none
    integer(ip), parameter :: n_targ = 4, n_max = 8
    real(dp) :: eig(n_max), evec(n, n_max), evec2(n, n_max), ref(n_targ), ax(n)
    real(dp), allocatable :: tmp(:)
    logical :: ok
    integer(ip) :: info, i, j
    integer :: n_failed = 0, it
    procedure(dgl_matvec), pointer :: metric_p
    real(dp) :: errs(8)
    logical :: oks(8)
!
    call setup()
    metric_p => mv_metric
!
! 1. Davidson, standard, random guess
    evec = 0; call dgl_davidson_driver(n, n_targ, n_max, mv_sym, pc_sym, eig, evec, ok, dgl_tol=1.0e-9_dp, dgl_info=info)
    call sorted(d, ref); call report("Davidson", ok .and. info == 0, eig, ref, residual_sym())
! 2. Davidson, generalized
    evec = 0; call dgl_davidson_driver(n, n_targ, n_max, mv_sym, pc_sym, eig, evec, ok, dgl_tol=1.0e-9_dp, &
                                       metvec=metric_p, dgl_info=info)
    call sorted(d/s, ref); call report("Davidson generalized", ok .and. info == 0, eig, ref, residual_gen())
! 3. LOBPCG, standard, simple guess, verbose with shift
    evec = 0; do i = 1, n_max; evec(i, i) = 1; end do
    call dgl_lobpcg_driver(n, n_targ, n_max, mv_sym, pc_spd, eig, evec, ok, dgl_tol=1.0e-9_dp, dgl_info=info, &
                           dgl_verbose=.false., dgl_shift=-100.0_dp, dgl_max_iter=300_ip)
    call sorted(d, ref); call report("LOBPCG", ok .and. info == 0, eig, ref, residual_sym())
! 4. LOBPCG generalized
    evec = 0; call dgl_lobpcg_driver(n, n_targ, n_max, mv_sym, pc_spd, eig, evec, ok, dgl_tol=1.0e-9_dp, &
                                     metvec=metric_p, dgl_info=info, dgl_max_iter=300_ip)
    call sorted(d/s, ref); call report("LOBPCG generalized", ok .and. info == 0, eig, ref, residual_gen())
! 5. non-symmetric, LR, biorthonormality
    evec = 0; evec2 = 0
    call dgl_davidson_nosym_driver(n, n_targ, n_max, mv_r, mv_l, pc_sym, "LR", eig, evec, ok, evec_2=evec2, &
                                   dgl_tol=1.0e-9_dp, dgl_info=info)
    call sorted(d, ref)
    call report("non-symmetric LR", ok .and. info == 0 .and. biortho_err() < 1.0e-10_dp, eig, ref, residual_nosym())
! 6. SMO-GD
    block
        real(dp) :: ev2(2*n, n_max)
        ev2 = 0
        call dgl_smogd_driver(2*n, n_targ, n_max, apb, amb, ident, ident, lrprec, eig, ev2, ok, dgl_tol=1.0e-9_dp, &
                              dgl_info=info, dgl_max_iter=200_ip)
        call sorted(sqrt(a*b), ref)
        call report("SMO-GD", ok .and. info == 0, eig, ref, 0.0_dp)
    end block
! 7. errors are reported, not stopping
    call dgl_davidson_driver(n, n_max + 1, n_max, mv_sym, pc_sym, eig, evec, ok, dgl_info=info)
    call check("error: n_targ > n_max", info == dgl_err_input .and. .not. ok)
    nullify (metric_p)
    call dgl_lobpcg_driver(n, n_targ, n_max, mv_sym, pc_spd, eig, evec, ok, metvec=metric_p, dgl_info=info)
    call check("error: disassociated metvec", info == dgl_err_input)
    call dgl_davidson_nosym_driver(n, n_targ, n_max, mv_r, mv_l, pc_sym, "LR", eig, evec, ok, dgl_info=info)
    call check("error: LR without evec_2", info == dgl_err_input)
    call dgl_davidson_driver(n, n_targ, n_max, zero_mv, pc_sym, eig, evec, ok, dgl_info=info, dgl_max_iter=5_ip)
    call check("zero matrix: no crash", .true.)
    write (*, "(a,i0)") "   (zero matrix: info = ", info
    call dgl_davidson_driver(n, n_targ, n_max, mv_sym, pc_sym, eig, evec, ok, dgl_memory=1_ip, dgl_memory_unit="KB", &
                             dgl_info=info)
    call check("error: memory", info == dgl_err_memory)
! 8. nested call
    evec = 0; call dgl_davidson_driver(n, n_targ, n_max, mv_nested, pc_sym, eig, evec, ok, dgl_tol=1.0e-9_dp, dgl_info=info)
    call sorted(d, ref); call report("Davidson, DiagLib called in matvec", ok .and. info == 0 .and. nested_calls > 0, &
                                      eig, ref, residual_sym())
! 9. concurrent calls from OpenMP threads, with different matrices
    !$omp parallel do num_threads(4) private(it, eig, evec, ok, info, ref) schedule(dynamic)
    do it = 1, 8
        scale = real(it, dp)
        evec = 0
        if (mod(it, 2) == 0) then
            call dgl_davidson_driver(n, n_targ, n_max, mv_nested, pc_sym, eig, evec, ok, dgl_tol=1.0e-9_dp, dgl_info=info)
        else
            call dgl_lobpcg_driver(n, n_targ, n_max, mv_sym, pc_spd, eig, evec, ok, dgl_tol=1.0e-9_dp, dgl_info=info, &
                                   dgl_max_iter=300_ip)
        end if
        call sorted(scale*d, ref)
        oks(it) = ok .and. info == 0
        errs(it) = maxval(abs(eig(:n_targ) - ref))
    end do
    !$omp end parallel do
    call check("8 concurrent calls (4 threads)", all(oks) .and. maxval(errs) < 1.0e-8_dp)
    write (*, "(a,es9.2)") "   max eigenvalue error of the concurrent calls: ", maxval(errs)
!
    write (*, "(a,i0,a)") "Fortran consumer: ", n_failed, " failures"
    if (n_failed > 0) error stop 1
contains
    subroutine sorted(x, r)
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: r(n_targ)
        real(dp) :: c(n)
        integer(ip) :: k, l
        c = x
        do k = 1, n_targ
            l = minloc(c, 1)
            r(k) = c(l)
            c(l) = huge(1.0_dp)
        end do
    end subroutine sorted
    real(dp) function residual_sym()
        integer(ip) :: k
        residual_sym = 0
        do k = 1, n_targ
            call mv_sym(n, 1_ip, evec(:, k), ax)
            residual_sym = max(residual_sym, norm2(ax - eig(k)*evec(:, k))/norm2(evec(:, k)))
        end do
    end function residual_sym
    real(dp) function residual_gen()
        integer(ip) :: k
        real(dp) :: bx(n)
        residual_gen = 0
        do k = 1, n_targ
            call mv_sym(n, 1_ip, evec(:, k), ax)
            call mv_metric(n, 1_ip, evec(:, k), bx)
            residual_gen = max(residual_gen, norm2(ax - eig(k)*bx)/norm2(evec(:, k)))
        end do
    end function residual_gen
    real(dp) function residual_nosym()
        integer(ip) :: k
        residual_nosym = 0
        do k = 1, n_targ
            call mv_r(n, 1_ip, evec(:, k), ax)
            residual_nosym = max(residual_nosym, norm2(ax - eig(k)*evec(:, k))/norm2(evec(:, k)))
            call mv_l(n, 1_ip, evec2(:, k), ax)
            residual_nosym = max(residual_nosym, norm2(ax - eig(k)*evec2(:, k))/norm2(evec2(:, k)))
        end do
    end function residual_nosym
    real(dp) function biortho_err()
        integer(ip) :: k, l
        biortho_err = 0
        do k = 1, n_targ
            do l = 1, n_targ
                biortho_err = max(biortho_err, abs(dot_product(evec2(:, k), evec(:, l)) - merge(1.0_dp, 0.0_dp, k == l)))
            end do
        end do
    end function biortho_err
    subroutine report(label, cond, e, r, res)
        character(len=*), intent(in) :: label
        logical, intent(in) :: cond
        real(dp), intent(in) :: e(:), r(:), res
        real(dp) :: err
        err = maxval(abs(e(:n_targ) - r))
        write (*, "(a40,a,es9.2,a,es9.2)") label, "  eig err ", err, "  residual ", res
        call check(label, cond .and. err < 1.0e-8_dp .and. res < 1.0e-6_dp)
    end subroutine report
    subroutine check(label, cond)
        character(len=*), intent(in) :: label
        logical, intent(in) :: cond
        if (.not. cond) n_failed = n_failed + 1
        write (*, "(a40,a)") label, merge("  PASSED", "  FAILED", cond)
    end subroutine check
end program ext_fortran
