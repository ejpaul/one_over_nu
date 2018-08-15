module grids_mod

	use stel_kinds
	use input_mod, only: nalpha, nlambda

	! Grids
	real(dp), dimension(:), allocatable :: alphas, lambdas
	! Spacing for alpha and lambda grids
	real(dp) :: dalpha
	real(dp), dimension(:), allocatable :: dlambda

	contains

! ===================================================
! subroutine init_grids
!
! This subroutine allocates and initializes the alphas
! and lambdas grids. 
!
! ===================================================
	subroutine init_grids

		use stel_constants
		use input_mod, only: nsurf

		implicit none

		integer :: ialpha, ilambda

		allocate(lambdas(nlambda))
		allocate(alphas(nalpha))
		allocate(dlambda(nsurf))

		do ialpha=1,nalpha
			alphas(ialpha) = twopi*(ialpha-1)/nalpha
		end do

		! Here lambda is scaled between 0 and 1
		do ilambda=1,nlambda
			lambdas(ilambda) = 1.0*(ilambda-1)/(nlambda-1)
		end do

		dlambda = lambdas(2)-lambdas(1)
		dalpha = alphas(2)-alphas(1)

	end subroutine init_grids

end module grids_mod
