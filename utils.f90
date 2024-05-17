module utils
  use real_precision
  implicit none
  real(dp),    allocatable :: a(:,:), s(:,:), apb(:,:), amb(:,:), aa(:,:), bb(:,:), sigma(:,:), delta(:,:), spd(:,:), smd(:,:)
!
! old complex implementation 
!
  complex(dp), allocatable :: caa(:,:), cbb(:,:), csigma(:,:), cdelta(:,:) !complex A, B, sigma, delta  matrices
  complex(dp), allocatable :: ca(:,:), cs(:,:) !complex Lambda, Omega matrice :: ca(:,:), cs(:,:) !complex Lambda, Omega matrices 
  complex(dp), allocatable :: capb(:,:), camb(:,:) !complex A+B, A-B matrices
! 
! new complex implementation 
!
  real(dp),    allocatable :: are(:,:), aim(:,:), sre(:,:), sim(:,:), apbre(:,:), apbim(:,:), ambre(:,:), ambim(:,:)
  real(dp),    allocatable :: aare(:,:), aaim(:,:), bbre(:,:), bbim(:,:), sigmare(:,:), sigmaim(:,:), deltare(:,:), deltaim(:,:)
  real(dp),    allocatable :: spdre(:,:), spdim(:,:), smdre(:,:), smdim(:,:)
  real(dp),    allocatable :: fourmat(:,:)
  real(dp), allocatable    :: lambdafull(:,:), omegafull(:,:)
!
  integer                  :: nmult, nsmult
  logical                  :: tdscf
  integer                  :: i_alg 
end module utils
!
