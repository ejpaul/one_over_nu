module splines_mod

	use ezspline_obj
	use ezspline
	use stel_kinds
	use input_mod, only: ntheta_spline, nzeta_spline
	use constants_mod

	real(dp), dimension(:), allocatable :: thetas_spline, zetas_spline
	real(dp), dimension(:,:,:), allocatable :: B_for_spline, dBdtheta_for_spline, dBdzeta_for_spline
	real(dp) :: dtheta_spline, dzeta_spline
	type(EZspline2_r8), dimension(:), allocatable :: spline_B, spline_dBdtheta, spline_dBdzeta

	contains

	! ===================================================
	! Subroutine init_spline
	!
	! This subroutine initializes the grids for the
	! cubic spline representation of B, dBdtheta, and
	! dBdzeta using the EZspline interface.
	!
	! ===================================================
	subroutine init_spline

		use geometry_mod, only: compute_geometry, nfp
		use input_mod, only: nsurf
		use constants_mod
		use stel_constants

		implicit none

		integer :: isurf, itheta, izeta, ierr
		real(dp), dimension(geometry_length) :: geometry
		real(dp) :: BB, theta, zeta, dBBdtheta, dBBdzeta

		allocate(B_for_spline(nsurf,ntheta_spline,nzeta_spline))
		allocate(dBdtheta_for_spline(nsurf,ntheta_spline,nzeta_spline))
		allocate(dBdzeta_for_spline(nsurf,ntheta_spline,nzeta_spline))
		allocate(thetas_spline(ntheta_spline))
		allocate(zetas_spline(nzeta_spline))

		! Construct grid for splines - must include periodic end points
		do itheta=1,ntheta_spline
			thetas_spline(itheta) = (2*pi)*(itheta-1)/(ntheta_spline-1)
		end do
		do izeta=1,nzeta_spline
			zetas_spline(izeta) = (2*pi/nfp)*(izeta-1)/(nzeta_spline-1)
		end do
		dtheta_spline = thetas_spline(2)-thetas_spline(1)
		dzeta_spline = zetas_spline(2)-zetas_spline(1)

		! Compute B on grid for splines
		do isurf=1,nsurf
			do itheta=1,ntheta_spline
				do izeta=1,nzeta_spline
					geometry = compute_geometry(isurf,thetas_spline(itheta),zetas_spline(izeta))
					B_for_spline(isurf,itheta,izeta) = geometry(B_index)
					dBdtheta_for_spline(isurf,itheta,izeta) = geometry(dBdtheta_index)
					dBdzeta_for_spline(isurf,itheta,izeta) = geometry(dBdzeta_index)
				end do
			end do
		end do

		allocate(spline_B(nsurf))
		allocate(spline_dBdtheta(nsurf))
		allocate(spline_dBdzeta(nsurf))
		do isurf=1,nsurf
			call ezspline_init(spline_B(isurf),ntheta_spline,nzeta_spline,(/-1,-1/),(/-1,-1/),ierr)
			call ezspline_error(ierr)
			spline_B(isurf)%x1 = thetas_spline
			spline_B(isurf)%x2 = zetas_spline

			call ezspline_setup(spline_B(isurf),B_for_spline(isurf,:,:),ierr,.true.)
			call ezspline_error(ierr)

			call ezspline_init(spline_dBdtheta(isurf),ntheta_spline,nzeta_spline,(/-1,-1/),(/-1,-1/),ierr)
			call ezspline_error(ierr)
			spline_dBdtheta(isurf)%x1 = thetas_spline
			spline_dBdtheta(isurf)%x2 = zetas_spline

			call ezspline_setup(spline_dBdtheta(isurf),dBdtheta_for_spline(isurf,:,:),ierr,.true.)
			call ezspline_error(ierr)

			call ezspline_init(spline_dBdzeta(isurf),ntheta_spline,nzeta_spline,(/-1,-1/),(/-1,-1/),ierr)
			call ezspline_error(ierr)
			spline_dBdzeta(isurf)%x1 = thetas_spline
			spline_dBdzeta(isurf)%x2 = zetas_spline

			call ezspline_setup(spline_dBdzeta(isurf),dBdzeta_for_spline(isurf,:,:),ierr,.true.)
			call ezspline_error(ierr)

		end do

		! Testing splines
		do isurf=1,nsurf
			do itheta=1,ntheta_spline
				do izeta=1,nzeta_spline
					theta = thetas_spline(itheta)+twopi
					zeta = zetas_spline(izeta)-twopi/nfp
					call spline_modulo(spline_B(isurf),theta,zeta)
					call ezspline_interp(spline_B(isurf),theta,zeta,BB,ierr)
					call ezspline_interp(spline_dBdtheta(isurf),theta,zeta,dBBdtheta,ierr)
					call ezspline_interp(spline_dBdzeta(isurf),theta,zeta,dBBdzeta,ierr)
					if (abs(BB-B_for_spline(isurf,itheta,izeta))>1.0e-6 .or. &
							abs(dBBdtheta-dBdtheta_for_spline(isurf,itheta,izeta))>1.0e-6 .or. &
							abs(dBBdzeta-dBdzeta_for_spline(isurf,itheta,izeta))>1.0e-6) then
						print *,"Error in spline is > 1e-6."
						stop
					end if
				end do
			end do
		end do

	end subroutine init_spline

! ===================================================
! Subroutine spline_modulo
!
! This subroutine should function as
! Ezspline_modulo, but allows negative arguments.
!
!
! Inputs:
! spline	ezspline object (type ezspline2_r8 assumed here)
!	theta		poloidal angle for evaluation
! zeta		toroidal angle for evaluation
!
! Outputs:
! theta		poloidal angle shifted to range of periodic
!					grid
! zeta		toroidal angle shifted to range of periodic
!					grid
! ===================================================
	subroutine spline_modulo(spline,theta,zeta)

		type(EZspline2_r8), intent(in) :: spline
		real(dp), intent(inout) :: theta, zeta

		! Check for periodicity
		if (spline%ibctype1(1)==-1 .and. spline%ibctype1(2)==-1) then
			do while (theta < spline%x1min)
				theta = theta + (spline%x1max-spline%x1min)
			end do
			do while (theta > spline%x1max)
				theta = theta - (spline%x1max-spline%x1min)
			end do
		end if
		! Check for periodicity
		if (spline%ibctype2(1)==-1 .and. spline%ibctype2(2)==-1) then
			do while (zeta < spline%x2min)
				zeta = zeta + (spline%x2max-spline%x2min)
			end do
			do while (zeta > spline%x2max)
				zeta = zeta - (spline%x2max-spline%x2min)
			end do
		end if
		! Check if within grid
		if (theta < spline%x1min .or. theta > spline%x1max) then
			print *,"Error in spline_modulo! "
		end if
		if (zeta < spline%x2min .or. zeta > spline%x2max) then
			print *,"Error in spline_modulo!"
		end if

	end subroutine spline_modulo

! ===================================================
! Function compute_B_spline
!
! This function computes the local magnetic field
! using the EZspline interface.
!
! Input:
! isurf		index in ssurf for calculation.
! theta		theta for calculation
! zeta		zeta for calculation.
!
! Output:
! compute_B_spline			interpolated B at specified location.
!
! ===================================================
	function compute_B_spline(isurf,theta,zeta)

		integer, intent(in) :: isurf
		real(dp) :: theta, zeta
		real(dp) :: compute_B_spline
		real(dp) :: zeta_eval, theta_eval
		integer :: ierr

		zeta_eval = zeta
		theta_eval = theta
		call spline_modulo(spline_B(isurf),theta_eval,zeta_eval)
		call ezspline_interp(spline_B(isurf),theta_eval,zeta_eval,compute_B_spline,ierr)
		call ezspline_error(ierr)

	end function compute_B_spline

! ===================================================
! Function compute_geometry_spline
!
! This function computes the local B, dBdtheta, and
! dBdzeta using the EzSpline interface.
!
! Input:
! isurf		index in ssurf for calculation.
! theta		theta for calculation
! zeta		zeta for calculation.
!
! Output:
! compute_geometry_spline		array of length geometry_length=4
!														This array can be indexed into
!														using B_index=1, dBdtheta=2,
!														dBdzeta=3 (defined in constants_mod).
!
! ===================================================
	function compute_geometry_spline(isurf,theta,zeta)

		integer, intent(in) :: isurf
		real(dp), intent(in) :: theta, zeta
		real(dp), dimension(geometry_length) :: compute_geometry_spline
		real(dp) :: zeta_eval, theta_eval
		integer :: ierr

		zeta_eval = zeta
		theta_eval = theta
		call spline_modulo(spline_B(isurf),theta_eval,zeta_eval)
		call ezspline_interp(spline_B(isurf),theta_eval,zeta_eval,compute_geometry_spline(B_index),ierr)
		call ezspline_error(ierr)
		call ezspline_interp(spline_dBdtheta(isurf),theta_eval,zeta_eval,&
				compute_geometry_spline(dBdtheta_index),ierr)
		call ezspline_error(ierr)
		call ezspline_interp(spline_dBdzeta(isurf),theta_eval,zeta_eval,&
				compute_geometry_spline(dBdzeta_index),ierr)
		call ezspline_error(ierr)

	end function compute_geometry_spline

end module splines_mod
