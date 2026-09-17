module dgl_interface
!* Public interface of DiagLib: kinds, error codes, interfaces of the user-supplied routines
! and of the drivers. This is the only module needed (and installed) to use DiagLib from
! Fortran: the drivers are implemented in submodules, and all the other modules are internal.
!
! DiagLib keeps no global state: the drivers can be called from inside the user-supplied
! routines of another driver call, and from different threads at the same time, as long as
! the user-supplied routines allow it.
    implicit none
!
#ifdef DGL_INT_KIND_8
    integer, parameter :: dgl_int = selected_int_kind(15)
!! Kind of the integers used by DiagLib (64-bit)
#elif DGL_INT_KIND_4
    integer, parameter :: dgl_int = selected_int_kind(8)
!! Kind of the integers used by DiagLib (32-bit)
#endif
    integer, parameter :: dgl_real = selected_real_kind(15)
!! Kind of the reals used by DiagLib (double precision)
!
! error codes returned in dgl_info
!
    integer(dgl_int), parameter :: dgl_success = 0_dgl_int
!! No error
    integer(dgl_int), parameter :: dgl_err_input = -1_dgl_int
!! Invalid input arguments
    integer(dgl_int), parameter :: dgl_err_memory = -2_dgl_int
!! Allocation failure or memory limit exceeded
    integer(dgl_int), parameter :: dgl_err_lapack = -3_dgl_int
!! A Lapack routine failed
    integer(dgl_int), parameter :: dgl_err_ortho = -4_dgl_int
!! An orthogonalization procedure failed
    integer(dgl_int), parameter :: dgl_err_mismatch = -5_dgl_int
!! Left and right eigenvalues of the non-symmetric driver do not match
!
! interfaces of the user-supplied routines
!
    abstract interface
        subroutine dgl_matvec(n, m, x, y)
!! Interface for the external routine performing the matrix-vector (or metric-vector) product
            import :: dgl_int, dgl_real
            implicit none
            integer(dgl_int), intent(in) :: n
!! Lenght of the vectors to multiply
            integer(dgl_int), intent(in) :: m
!! Number of vectors to multiply
            real(dgl_real), dimension(n, m), intent(in) :: x
!! Input vectors
            real(dgl_real), dimension(n, m), intent(inout) :: y
!! Output vectors
        end subroutine dgl_matvec
!
        subroutine dgl_precnd(n, m, shift, x, y)
!! Interface for the external routine performing the preconditioning of the residuals
            import :: dgl_int, dgl_real
            implicit none
            integer(dgl_int), intent(in) :: n
!! Lenght of the vectors
            integer(dgl_int), intent(in) :: m
!! Number of vectors
            real(dgl_real), intent(in) :: shift
!! Shift for the preconditioner: minus the lowest non-converged approximate eigenvalue,
!! or zero (see the dgl_precnd_shift argument of the drivers)
            real(dgl_real), dimension(n, m), intent(in) :: x
!! Input vectors
            real(dgl_real), dimension(n, m), intent(inout) :: y
!! Output vectors
        end subroutine dgl_precnd
!
        subroutine dgl_smogd_precnd(n, m, fac, xp, xm, yp, ym)
!! Interface for the external routine performing the preconditioning in SMO-GD
            import :: dgl_int, dgl_real
            implicit none
            integer(dgl_int), intent(in) :: n
!! Lenght of the vectors
            integer(dgl_int), intent(in) :: m
!! Number of vectors
            real(dgl_real), intent(in) :: fac
!! Inverse of the current approximation to the eigenvalue
            real(dgl_real), dimension(n, m), intent(in) :: xp
!! Input vectors, plus combination
            real(dgl_real), dimension(n, m), intent(in) :: xm
!! Input vectors, minus combination
            real(dgl_real), dimension(n, m), intent(inout) :: yp
!! Output vectors, plus combination
            real(dgl_real), dimension(n, m), intent(inout) :: ym
!! Output vectors, minus combination
        end subroutine dgl_smogd_precnd
    end interface
!
! drivers
!
    interface
        module subroutine dgl_davidson_driver(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                                   dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                   dgl_shift, dgl_memory, dgl_memory_unit, metvec, dgl_info, dgl_precnd_shift)
!! # Driver for Davidson-Liu symmetric diagonalization
!! Can solve both standard and generalized eigenvalue problems.
!! In the latter case you need to pass the optional argument [[metvec]] as a pointer to your routine.
!! @note
!! eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
!! @endnote
            implicit none
            integer(dgl_int), intent(in) :: n
!! Size of the matrix to be diagonalized
            integer(dgl_int), intent(in) :: n_targ
!! Number of required eigenpairs.
            integer(dgl_int), intent(in) :: n_max
!! Maximum size of the search space. Should be >= n_targ.
            real(dgl_real), dimension(n_max), intent(inout) :: eig
!! Computed eigenvalues
            real(dgl_real), dimension(n, n_max), intent(inout) :: evec
!! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
            logical, intent(inout) :: ok
!! True if davidson converged
            procedure(dgl_matvec) :: matvec
!! External subroutine that performs the matrix-vector multiplication
            procedure(dgl_precnd) :: precnd
!! External subroutine that applies a preconditioner
            logical, optional, intent(in) :: dgl_verbose
!! Verbose mode. Default = .false.
            integer(dgl_int), optional, intent(in) :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
            integer(dgl_int), optional, intent(in) :: dgl_dav_iter
!! Maximum number of iterations before Davidson restart. Default = \(25\)
            integer(dgl_int), optional, intent(in) :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
            character(len=2), optional, intent(in) :: dgl_memory_unit
!! Unit of memory. Default = MB
            real(dgl_real), optional, intent(in) :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
            real(dgl_real), optional, intent(in) :: dgl_shift
!! Constant added to the eigenvalues when they are printed, e.g., the energy of the frozen core
!! to print total energies. It does not affect the computation nor the returned eigenvalues.
!! Default = \(0.\)
            procedure(dgl_matvec), pointer, optional :: metvec
!! Pointer to External subroutine that applies the metric-vector multiplication
            integer(dgl_int), optional, intent(out) :: dgl_info
!! Error status: dgl_success (0) or one of the (negative) dgl_err_* codes.
!! If not present, DiagLib stops the program when an error occurs.
            logical, optional, intent(in) :: dgl_precnd_shift
!! If true, the shift passed to precnd is minus the lowest non-converged eigenvalue, as in
!! Davidson's method; if false, the shift is zero. Default = .true.
        end subroutine dgl_davidson_driver

        module subroutine dgl_lobpcg_driver(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                                 dgl_verbose, dgl_max_iter, dgl_tol, &
                                 dgl_shift, dgl_memory, dgl_memory_unit, metvec, dgl_info, dgl_precnd_shift)
!! # Driver for LOBPCG symmetric diagonalization
!! Can solve both standard and generalized eigenvalue problems.
!! In the latter case you need to pass the optional argument [[metvec]] as a pointer to your routine.
!! @note
!! eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
!! @endnote
            implicit none
            integer(dgl_int), intent(in) :: n
!! Size of the matrix to be diagonalized
            integer(dgl_int), intent(in) :: n_targ
!! Number of required eigenpairs.
            integer(dgl_int), intent(in) :: n_max
!! Maximum size of the search space. Should be >= n_targ.
            real(dgl_real), dimension(n_max), intent(inout) :: eig
!! Computed eigenvalues
            real(dgl_real), dimension(n, n_max), intent(inout) :: evec
!! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
            logical, intent(inout) :: ok
!! True if davidson converged
            procedure(dgl_matvec) :: matvec
!! External subroutine that performs the matrix-vector multiplication
            procedure(dgl_precnd) :: precnd
!! External subroutine that applies a preconditioner
            logical, optional, intent(in) :: dgl_verbose
!! Verbose mode. Default = .false.
            integer(dgl_int), optional, intent(in) :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
            integer(dgl_int), optional, intent(in) :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
            character(len=2), optional, intent(in) :: dgl_memory_unit
!! Unit of memory. Default = MB
            real(dgl_real), optional, intent(in) :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
            real(dgl_real), optional, intent(in) :: dgl_shift
!! Constant added to the eigenvalues when they are printed, e.g., the energy of the frozen core
!! to print total energies. It does not affect the computation nor the returned eigenvalues.
!! Default = \(0.\)
            procedure(dgl_matvec), pointer, optional :: metvec
!! Pointer to External subroutine that applies the metric-vector multiplication
            integer(dgl_int), optional, intent(out) :: dgl_info
!! Error status: dgl_success (0) or one of the (negative) dgl_err_* codes.
!! If not present, DiagLib stops the program when an error occurs.
            logical, optional, intent(in) :: dgl_precnd_shift
!! If true, the shift passed to precnd is minus the lowest non-converged eigenvalue; if false,
!! the shift is zero. LOBPCG works best with a positive definite preconditioner that stays well
!! conditioned (e.g., an approximation of the inverse of the matrix), which a preconditioner
!! shifted by an approximate eigenvalue is not. Default = .false.
        end subroutine dgl_lobpcg_driver

        module subroutine dgl_davidson_nosym_driver(n, n_targ, n_max, matvec_r, matvec_l, precnd, side, &
                                         eig, evec_1, ok, evec_2, &
                                         dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                         dgl_shift, dgl_memory, dgl_memory_unit, dgl_info, dgl_precnd_shift)
!! # Driver for Davidson-Liu non-symmetric diagonalization
!! Non-symmetric davidson diagonalization is commonly encountered in EOM-CC theory.
!! This driver can eveluate both Left and Right eigenvectors.
!! Only standard eigenvalue problems.
!! @note
!! eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
!! @endnote
!!
!!
!! Lapack does not order the eigenvalues of the reduced (non-symmetric) matrix: the ritz pairs
!! are sorted by increasing real part, and the n_max lowest real ones are kept.
!! To follow each root through the iterations, also in case of root flipping or near
!! degeneracies, the kept ritz pairs are then assigned to the slots (the positions in eig and
!! the columns of the eigenvectors) of the previous iteration by maximum overlap.
!! Convergence is checked on the n_targ lowest roots, whatever slot they are in, and the
!! results are returned sorted by increasing eigenvalue.
            implicit none
            integer(dgl_int), intent(in) :: n
!! Size of the matrix to be diagonalized
            integer(dgl_int), intent(in) :: n_targ
!! Number of required eigenpairs.
            integer(dgl_int), intent(in) :: n_max
!! Maximum size of the search space. Should be >= n_targ
            character(len=2), intent(in) :: side
!! String to decide which eigenvectors to compute and whether to compute
!! both. Possible values are "R ", "L " or "LR"
            real(dgl_real), dimension(n_max), intent(inout) :: eig
!! Computed eigenvalues
            real(dgl_real), dimension(n, n_max), intent(inout) :: evec_1
!! First set of computed vectors. In input it should contain a guess
!! If side="LR" contains the Right ones
            real(dgl_real), dimension(n, n_max), optional, intent(inout) :: evec_2
!! Second set of computed vectors. First set of converged vectors is used as guess.
!! If side="LR" contains the Left ones
!! If side="LR" and the driver converged, the first n_targ left and right eigenvectors are
!! biorthonormal: evec_2(:, i) . evec_1(:, j) = delta_ij, with normalized right eigenvectors.
            logical, intent(inout) :: ok
!! True if davidson converged
            procedure(dgl_matvec) :: matvec_r
!! External subroutine that performs the matrix-vector multiplication for
!! right eigenvectors
            procedure(dgl_matvec) :: matvec_l
!! External subroutine that performs the matrix-vector multiplication for
!! left eigenvectors
            procedure(dgl_precnd) :: precnd
!! External subroutine that applies a preconditioner
            logical, optional, intent(in) :: dgl_verbose
!! Verbose mode. Default = .false.
            integer(dgl_int), optional, intent(in) :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
            integer(dgl_int), optional, intent(in) :: dgl_dav_iter
!! Maximum number of iterations before Davidson restart. Default = \(25\)
            integer(dgl_int), optional, intent(in) :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
            character(len=2), optional, intent(in) :: dgl_memory_unit
!! Unit of memory. Default = MB
            real(dgl_real), optional, intent(in) :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
            real(dgl_real), optional, intent(in) :: dgl_shift
!! Constant added to the eigenvalues when they are printed, e.g., the energy of the frozen core
!! to print total energies. It does not affect the computation nor the returned eigenvalues.
!! Default = \(0.\)
            integer(dgl_int), optional, intent(out) :: dgl_info
!! Error status: dgl_success (0) or one of the (negative) dgl_err_* codes.
!! If not present, DiagLib stops the program when an error occurs.
            logical, optional, intent(in) :: dgl_precnd_shift
!! If true, the shift passed to precnd is minus the lowest non-converged eigenvalue, as in
!! Davidson's method; if false, the shift is zero. Default = .true.
        end subroutine dgl_davidson_nosym_driver

        module subroutine dgl_smogd_driver(n2, n_targ, n_max, apbmul, ambmul, &
                                spdmul, smdmul, lrprec, eig, evec, ok, &
                                dgl_verbose, dgl_tol, dgl_max_iter, dgl_dav_iter, &
                                dgl_memory, dgl_memory_unit, dgl_info)
!!# Driver for the efficient solution to the Linear-Response CASSCF problem
!! \begin{equation}
!! \begin{bmatrix} \begin{pmatrix}
!! A & B \\
!! B & A
!! \end{pmatrix}
!! -
!! \omega
!! \begin{pmatrix}
!! S & D \\
!! -D & -S
!! \end{pmatrix} \end{bmatrix}
!! \begin{pmatrix}
!! Y \\
!! Z
!! \end{pmatrix}
!! =
!! \begin{pmatrix}
!! 0 \\
!! 0
!! \end{pmatrix},
!! \label{eq:respeq}
!! \end{equation}
!!
!! Where A, B, S are symmetric matrices and D is antysimmetric.
!!
!! If \(\begin{bmatrix} w, \begin{pmatrix} Y \\ Z \end{pmatrix} \end{bmatrix}\) are a solution,
!! then \(\begin{bmatrix} -w, \begin{pmatrix} Z \\ Y \end{pmatrix} \end{bmatrix}\) is also a solution.
!!
!! Following J. Chem. Phys., 118, 522 (2003), we enforce this property in
!! the iterative procedure by expanding the eigenvector as
!!
!! \begin{equation}
!! \begin{pmatrix} Y \\ Z \end{pmatrix} =
!! \begin{pmatrix} b^+ \\ b^+ \end{pmatrix} +
!! \begin{pmatrix} b^- \\ -b^- \end{pmatrix}
!! \end{equation}
!!
!! This routine performs the Swapped Metric-Orthogonal -- Generalized Davidsion,
!! therefore solves the associate problem:
!!
!!\begin{equation}
!! \begin{bmatrix}
!! \begin{pmatrix}
!! S & D \\
!! -D & -S
!! \end{pmatrix}
!! -
!! \frac{1}{\omega}
!! \begin{pmatrix}
!! A & B \\
!! B & A
!! \end{pmatrix}
!! \end{bmatrix}
!! \begin{pmatrix}
!! Y \\
!! Z
!! \end{pmatrix}
!! =
!! \begin{pmatrix}
!! 0 \\
!! 0
!! \end{pmatrix},
!! \label{eq:respeq_smogd}
!!\end{equation}
!!
!! using the casida matrix, which is symmetric and positive definite, as
!! the metric. This allows us to use expansion vectors that are orthogonal
!! with respect to the dot product defined by the metric, which in turn
!! results in a Rayleigh-Ritz procedure that requires the solution of a
!! symmetric standard eigenvalue problem
!!
!! \begin{equation}
!! \begin{pmatrix}
!! 0 & s^T \\
!! s & 0
!! \end{pmatrix}
!! \begin{pmatrix}
!! u^+ \\
!! u^-
!! \end{pmatrix}
!! =
!! \frac{1}{\omega}
!! \begin{pmatrix}
!! u^+ \\
!! u^-
!! \end{pmatrix},
!! \label{eq:krylov}
!! \end{equation}
!!
!! which can be reduced to a half-sized eigenvalue problem
!!
!! \(s^T s u^+ = \left(\frac{1}{\omega}\right)^2 u^+ \\
!! u^- = \frac{1}{\omega} Su^+\)
!!
!! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
            implicit none
            integer(dgl_int), intent(in) :: n2
!! Totals size of the generalized eigenvalue
            integer(dgl_int), intent(in) :: n_targ
!! Number of required eigenpairs.
            integer(dgl_int), intent(in) :: n_max
!! Maximum size of the search space. Should be >= n_targ.
            real(dgl_real), dimension(n_max), intent(inout) :: eig
!! Computed eigenvalues
            real(dgl_real), dimension(n2, n_max), intent(inout) :: evec
!! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
            logical, intent(inout) :: ok
!! True if davidson converged
            procedure(dgl_matvec) :: apbmul
!! External subroutine that performs the matrix-vector multiplication with A+B
            procedure(dgl_matvec) :: ambmul
!! External subroutine that performs the matrix-vector multiplication with A-B
            procedure(dgl_matvec) :: spdmul
!! External subroutine that performs the matrix-vector multiplication with S+D
            procedure(dgl_matvec) :: smdmul
!! External subroutine that performs the matrix-vector multiplication with S-D
            procedure(dgl_smogd_precnd) :: lrprec
!! External subroutine that applies a preconditioner to both plus and minus vectors
            logical, optional, intent(in) :: dgl_verbose
!! Verbose mode. Default = .false.
            integer(dgl_int), optional, intent(in) :: dgl_dav_iter
!! Maximum number of iterations before Davidson restart. Default = \(25\)
            integer(dgl_int), optional, intent(in) :: dgl_max_iter
!! Maximum number of allowed iterations. Default = \(100\)
            integer(dgl_int), optional, intent(in) :: dgl_memory
!! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
            character(len=2), optional, intent(in) :: dgl_memory_unit
!! Unit of memory. Default = MBs
            real(dgl_real), optional, intent(in) :: dgl_tol
!! Convergence threshold on residuals norms. Default = \(10^{-7}\)
            integer(dgl_int), optional, intent(out) :: dgl_info
!! Error status: dgl_success (0) or one of the (negative) dgl_err_* codes.
!! If not present, DiagLib stops the program when an error occurs.
        end subroutine dgl_smogd_driver
    end interface
!
end module dgl_interface
