module dgl_drivers_interfaces
!* Module containing the explicit interfaces to the main drivers in DiagLib.
! This has to be included if one is interested in using the optional arguments
! (e.g. generalized problems).
    use dgl_global_utils, only: dp
    use dgl_external_interfaces
    implicit none

    interface
    !! Interface for main DiagLib drivers
        subroutine davidson_driver(n,n_targ,n_max,matvec,precnd,eig,evec,ok, &
              dgl_verbose, dgl_max_iter, dgl_tol, dgl_dav_iter, &
              dgl_shift, dgl_memory, dgl_memory_unit, metvec)
        !! ### Interface for Davidson-Liu symmetric diagonalization
        !! Can solve both standard and generalized eigenvalue problems.
        !! In the latter case you need to pass the optional argument [[metvec]] as a pointer to your routine.
        !! Moreover, you need to use the [[dgl_drivers_interfaces]] module included in this library.
        !!
        !! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
            import
            implicit none
            integer,                      intent(in)    :: n
            !! Size of the matrix to be diagonalized
            integer,                      intent(in)    :: n_targ
            !! Number of required eigenpairs.
            integer,                      intent(in)    :: n_max
            !! Maximum size of the search space. Should be >= n_targ.
            real(dp), dimension(n_max),   intent(inout) :: eig
            !! Computed eigenvalues    
            real(dp), dimension(n,n_max), intent(inout) :: evec
            !! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
            logical,                      intent(inout) :: ok
            !! True if davidson converged
            procedure(matvec_) :: matvec
            !! External subroutine that performs the matrix-vector multiplication
            procedure(precnd_) :: precnd
            !! External subroutine that applies a preconditioner
            logical,  optional,            intent(in)    :: dgl_verbose
            !! Verbose mode. Default = .false.
            integer,  optional,            intent(in)    :: dgl_max_iter
            !! Maximum number of allowed iterations. Default = \(100\)
            integer,  optional,            intent(in)    :: dgl_dav_iter
            !! Maximum number of iterations before Davidson restart. Default = \(25\)
            integer,  optional,            intent(in)    :: dgl_memory
            !! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
            character(len=2),  optional,   intent(in)    :: dgl_memory_unit
            !! Unit of memory. Default = MBs
            real(dp), optional,            intent(in)    :: dgl_tol
            !! Convergence threshold on residuals norms. Default = \(10^{-7}\)
            real(dp), optional,            intent(in)    :: dgl_shift
            !! Diagonal level shifting parameter. Default = \(0.\)
            procedure(metvec_), pointer, optional :: metvec
            !! Pointer to External subroutine that applies the metric-vector multiplication
        end subroutine davidson_driver

        subroutine lobpcg_driver(n,n_targ,n_max,matvec,precnd,eig,evec,ok, &
              dgl_verbose, dgl_max_iter, dgl_tol, &
              dgl_shift, dgl_memory,dgl_memory_unit, metvec)
        !! ### Interface for LOBPCG symmetric diagonalization
        !! Can solve both standard and generalized eigenvalue problems.
        !! In the latter case you need to pass the optional argument [[metvec]] as a pointer to your routine.
        !! Moreover, you need to use the [[dgl_drivers_interfaces]] module included in this library.
        !!
        !! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
            import
            implicit none
            integer,                      intent(in)    :: n
            !! Size of the matrix to be diagonalized
            integer,                      intent(in)    :: n_targ
            !! Number of required eigenpairs.
            integer,                      intent(in)    :: n_max
            !! Maximum size of the search space. Should be >= n_targ.
            real(dp), dimension(n_max),   intent(inout) :: eig
            !! Computed eigenvalues    
            real(dp), dimension(n,n_max), intent(inout) :: evec
            !! Computed eigenvectors. In input, it should contain a guess for the eigenvectors
            logical,                      intent(inout) :: ok
            !! True if davidson converged
            procedure(matvec_) :: matvec
            !! External subroutine that performs the matrix-vector multiplication
            procedure(precnd_) :: precnd
            !! External subroutine that applies a preconditioner
            logical,  optional,            intent(in)    :: dgl_verbose
            !! Verbose mode. Default = .false.
            integer,  optional,            intent(in)    :: dgl_max_iter
            !! Maximum number of allowed iterations. Default = \(100\)
            integer,  optional,            intent(in)    :: dgl_memory
            !! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
            character(len=2),  optional,   intent(in)    :: dgl_memory_unit
            !! Unit of memory. Default = MBs
            real(dp), optional,            intent(in)    :: dgl_tol
            !! Convergence threshold on residuals norms. Default = \(10^{-7}\)
            real(dp), optional,            intent(in)    :: dgl_shift
            !! Diagonal level shifting parameter. Default = \(0.\)
            procedure(metvec_), pointer, optional :: metvec
            !! Pointer to External subroutine that applies the metric-vector multiplication
        end subroutine lobpcg_driver
        
        subroutine davidson_nosym_driver(n,n_targ,n_max,matvec_r,matvec_l,precnd,side, &
              eig,evec_r,evec_l,ok,&
              dgl_verbose,dgl_tol,dgl_max_iter,dgl_dav_iter,&
              dgl_shift,dgl_memory,dgl_memory_unit)
        !! ### Interface for Davidson-Liu non-symmetric diagonalization
        !! Can solve solve only standard eigenvalue problems. Can eveluate both Left and Right eigenvectors.
        !! To pass any optional argument, you need to use [[dgl_drivers_interfaces]], a module included in this library.
        !!
        !! **Note:** eig and evec should be allocated (n_max) and (n,n_max), where \(n_{max} \ge n_{act}\).
              import
              implicit none
                integer,                      intent(in)    :: n
            !! Size of the matrix to be diagonalized
                integer,                      intent(in)    :: n_targ
            !! Number of required eigenpairs.
                integer,                      intent(in)    :: n_max
            !! Maximum size of the search space. Should be >= n_targ
                integer,                      intent(in)    :: side
            !! Integer to decide which eigenvectors to compute and whether to compute
            !! them togheter or separately
                real(dp), dimension(n_max),   intent(inout) :: eig
            !! Computed eigenvalues
                real(dp), dimension(n,n_max), intent(inout) :: evec_l
            !! Computed Left eigenvectors. In input, it should contain their guess
                real(dp), dimension(n,n_max), intent(inout) :: evec_r
            !! Computed Right eigenvectors. In input, it should contain their guess
                logical,                      intent(inout) :: ok
            !! True if davidson converged
                procedure(matvec_) :: matvec_r
            !! External subroutine that performs the matrix-vector multiplication for
            !! right eigenvectors
                procedure(matvec_) :: matvec_l
            !! External subroutine that performs the matrix-vector multiplication for
            !! left eigenvectors
                procedure(precnd_) :: precnd
            !! External subroutine that applies a preconditioner
                logical,  optional,            intent(in)    :: dgl_verbose
            !! Verbose mode. Default = .false.
                integer,  optional,            intent(in)    :: dgl_max_iter
            !! Maximum number of allowed iterations. Default = \(100\)
                integer,  optional,            intent(in)    :: dgl_dav_iter
            !! Maximum number of iterations before Davidson restart. Default = \(25\)
                integer,  optional,            intent(in)    :: dgl_memory
            !! Maximum memory that DiagLib is allowed to use. Default = \(80\)MBs
                character(len=2),  optional,   intent(in)    :: dgl_memory_unit
            !! Unit of memory. Default = MBs
                real(dp), optional,            intent(in)    :: dgl_tol
            !! Convergence threshold on residuals norms. Default = \(10^{-7}\)
                real(dp), optional,            intent(in)    :: dgl_shift
            !! Diagonal level shifting parameter. Default = \(0.\)
        end subroutine davidson_nosym_driver

    end interface

end module dgl_drivers_interfaces