module utils
  use real_precision
  implicit none
  real(dp), allocatable :: a(:,:), s(:,:), apb(:,:), amb(:,:), aa(:,:), &
    bb(:,:), sigma(:,:), delta(:,:), spd(:,:), smd(:,:), a_t(:,:)
  integer               :: nmult, nsmult
  logical               :: tdscf
  integer               :: i_alg 
end module utils
!
