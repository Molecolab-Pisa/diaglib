module utils
  use real_precision
  implicit none
  real(dp), allocatable :: a(:,:), s(:,:), apb(:,:), amb(:,:), aa(:,:), bb(:,:), sigma(:,:), delta(:,:), spd(:,:), smd(:,:)
  integer               :: nmult, nsmult
end module utils
!
