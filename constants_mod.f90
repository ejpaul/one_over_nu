module constants_mod

	use stel_kinds

	! Fixed lengths and indices
	integer, parameter :: string_length = 100
	integer, parameter :: integrand_length = 6
	integer, parameter :: dKdalpha_index = 1
	integer, parameter :: H_index = 2
	integer, parameter :: I_index = 3
	integer, parameter :: J_index = 4
	integer, parameter :: d2Kdalpha2_index = 5
	integer, parameter :: dIdalpha_index = 6
	integer, parameter :: s_wish_length = 100

	! Fixed lengths and indices
	integer, parameter :: geometry_length = 5
	integer, parameter :: B_index = 1
	integer, parameter :: dBdtheta_index = 2
	integer, parameter :: Bdotgradzeta_index = 3
	integer, parameter :: dBdotgradzetadtheta_index = 4
	integer, parameter :: d2Bdotgradzetadtheta2_index = 5

	! Physical constants
	! permitivity of free space (F/m)
	real(dp) :: epsilon0 = 8.8542e-12
	! elementary charge (C)
	real(dp) :: e = 1.6022e-19

end module constants_mod
