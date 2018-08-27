module splines_mod

	use ezspline_obj
	use ezspline
	use stel_kinds
	use input_mod, only: ntheta_spline, nzeta_spline, geometry_option, output_P_tensor
	use constants_mod

	real(dp), dimension(:), allocatable :: thetas_spline, zetas_spline
	real(dp), dimension(:,:,:), allocatable :: B_for_spline, dBdtheta_for_spline, dBdzeta_for_spline, Bdotgradzeta_for_spline, dBdotgradzetadtheta_for_spline, d2Bdtheta2_for_spline, &
		d2Bdotgradzetadtheta2_for_spline
	real(dp) :: dtheta_spline, dzeta_spline
	type(EZspline2_r8), dimension(:), allocatable :: spline_B, spline_dBdtheta, spline_dBdzeta, &
	spline_d2Bdtheta2, spline_Bdotgradzeta, spline_dBdotgradzetadtheta, spline_d2Bdotgradzetadtheta2

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

		use geometry_mod, only: compute_geometry, nfp, compute_d2Bdtheta2, compute_dBdzeta, &
			compute_d2Bdotgradzetadtheta2
		use input_mod, only: nsurf, output_P_tensor, geometry_option
		use constants_mod
		use stel_constants

		implicit none

		integer :: isurf, itheta, izeta, ierr
		real(dp), dimension(geometry_length) :: geometry
		real(dp) :: BB, theta, zeta, dBBdtheta, dBBdzeta, Bdotgradzeta, dBdotgradzetadtheta

		allocate(B_for_spline(nsurf,ntheta_spline,nzeta_spline))
		allocate(dBdtheta_for_spline(nsurf,ntheta_spline,nzeta_spline))
		allocate(dBdzeta_for_spline(nsurf,ntheta_spline,nzeta_spline))
		allocate(thetas_spline(ntheta_spline))
		allocate(zetas_spline(nzeta_spline))
		if (output_P_tensor) then
			allocate(d2Bdtheta2_for_spline(nsurf,ntheta_spline,nzeta_spline))
			if (geometry_option==2) then
				allocate(d2Bdotgradzetadtheta2_for_spline(nsurf,ntheta_spline,nzeta_spline))
			end if
		end if
		if (geometry_option==2) then
			allocate(Bdotgradzeta_for_spline(nsurf,ntheta_spline,nzeta_spline))
			allocate(dBdotgradzetadtheta_for_spline(nsurf,ntheta_spline,nzeta_spline))
		end if

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
					dBdzeta_for_spline(isurf,itheta,izeta) = &
						compute_dBdzeta(isurf,thetas_spline(itheta),zetas_spline(izeta))
					if (output_P_tensor) then
						d2Bdtheta2_for_spline(isurf,itheta,izeta) = &
							compute_d2Bdtheta2(isurf,thetas_spline(itheta),zetas_spline(izeta))
						if (geometry_option == 2) then
							d2Bdotgradzetadtheta2_for_spline(isurf,itheta,izeta) = &
								compute_d2Bdotgradzetadtheta2(isurf,thetas_spline(itheta),zetas_spline(izeta))
						end if
					end if
					if (geometry_option==2) then
						Bdotgradzeta_for_spline(isurf,itheta,izeta) = geometry(Bdotgradzeta_index)
						dBdotgradzetadtheta_for_spline(isurf,itheta,izeta) = &
							geometry(dBdotgradzetadtheta_index)
					end if
				end do
			end do
		end do

		allocate(spline_B(nsurf))
		allocate(spline_dBdtheta(nsurf))
		allocate(spline_dBdzeta(nsurf))

		if (output_P_tensor) then
			allocate(spline_d2Bdtheta2(nsurf))
			if (geometry_option==2) then
				allocate(spline_d2Bdotgradzetadtheta2(nsurf))
			end if
		end if
		if (geometry_option==2) then
			allocate(spline_Bdotgradzeta(nsurf))
			allocate(spline_dBdotgradzetadtheta(nsurf))
		end if
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

			if (geometry_option == 2) then ! PEST - eventually I should change this to use transform grid
				call ezspline_init(spline_Bdotgradzeta(isurf),ntheta_spline,nzeta_spline,(/-1,-1/),&
					(/-1,-1/),ierr)
				call ezspline_error(ierr)
				spline_Bdotgradzeta(isurf)%x1 = thetas_spline
				spline_Bdotgradzeta(isurf)%x2 = zetas_spline
				call ezspline_setup(spline_Bdotgradzeta(isurf),Bdotgradzeta_for_spline(isurf,:,:),ierr,.true.)
				call ezspline_error(ierr)

				call ezspline_init(spline_dBdotgradzetadtheta(isurf),ntheta_spline,nzeta_spline,(/-1,-1/),&
					(/-1,-1/),ierr)
				call ezspline_error(ierr)
				spline_dBdotgradzetadtheta(isurf)%x1 = thetas_spline
				spline_dBdotgradzetadtheta(isurf)%x2 = zetas_spline
				call ezspline_setup(spline_dBdotgradzetadtheta(isurf),dBdotgradzetadtheta_for_spline(isurf,:,:),ierr,.true.)
				call ezspline_error(ierr)
			end if

			if (output_P_tensor) then
				call ezspline_init(spline_d2Bdtheta2(isurf),ntheta_spline,&
					nzeta_spline,(/-1,-1/),(/-1,-1/),ierr)
				call ezspline_error(ierr)
				spline_d2Bdtheta2(isurf)%x1 = thetas_spline
				spline_d2Bdtheta2(isurf)%x2 = zetas_spline
				call ezspline_setup(spline_d2Bdtheta2(isurf),&
					d2Bdtheta2_for_spline(isurf,:,:),ierr,.true.)
				call ezspline_error(ierr)
				if (geometry_option==2) then
					call ezspline_init(spline_d2Bdotgradzetadtheta2(isurf),ntheta_spline,&
						nzeta_spline,(/-1,-1/),(/-1,-1/),ierr)
					call ezspline_error(ierr)
					spline_d2Bdotgradzetadtheta2(isurf)%x1 = thetas_spline
					spline_d2Bdotgradzetadtheta2(isurf)%x2 = zetas_spline
					call ezspline_setup(spline_d2Bdotgradzetadtheta2(isurf),&
						d2Bdotgradzetadtheta2_for_spline(isurf,:,:),ierr,.true.)
					call ezspline_error(ierr)
				end if
			end if
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
					if (geometry_option==2) then
						call ezspline_interp(spline_Bdotgradzeta(isurf),theta,zeta,Bdotgradzeta,ierr)
						call ezspline_interp(spline_dBdotgradzetadtheta(isurf),theta,zeta,dBdotgradzetadtheta,ierr)
						if (abs(Bdotgradzeta-Bdotgradzeta_for_spline(isurf,itheta,izeta))>1.0e-6 .or. &
							abs(dBdotgradzetadtheta-dBdotgradzetadtheta_for_spline(isurf,itheta,izeta))>1.0e-6) then
							print *,"Error in spline is > 1e-6"
							stop
						end if
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
		call ezspline_interp(spline_B(isurf),theta_eval,zeta_eval,&
			compute_B_spline,ierr)
		call ezspline_error(ierr)

	end function compute_B_spline

! ===================================================
! Function compute_dBzeta_spline
!
! This function computes the local dBdzeta
! using the EZspline interface.
!
! Input:
! isurf		index in ssurf for calculation.
! theta		theta for calculation
! zeta		zeta for calculation.
!
! Output:
! compute_dBdzeta_spline	interpolated dBdzeta at specified location.
!
! ===================================================
	function compute_dBdzeta_spline(isurf,theta,zeta)

		integer, intent(in) :: isurf
		real(dp) :: theta, zeta
		real(dp) :: compute_dBdzeta_spline
		real(dp) :: zeta_eval, theta_eval
		integer :: ierr

		zeta_eval = zeta
		theta_eval = theta
		call spline_modulo(spline_dBdzeta(isurf),theta_eval,zeta_eval)
		call ezspline_interp(spline_dBdzeta(isurf),theta_eval,zeta_eval,&
			compute_dBdzeta_spline,ierr)
		call ezspline_error(ierr)

	end function compute_dBdzeta_spline

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
! compute_geometry_spline		array of length geometry_length
!														This array can be indexed into
!														using parameters defined in
!														constants_mod.
!
! ===================================================
	function compute_geometry_spline(isurf,theta,zeta)

		use geometry_mod, only: Boozer_G, Boozer_I, iota
		use stel_constants

		integer, intent(in) :: isurf
		real(dp), intent(in) :: theta, zeta
		real(dp), dimension(geometry_length) :: compute_geometry_spline
		real(dp) :: zeta_eval, theta_eval, BB, dBBdtheta, sign_Bdotgradzeta
		integer :: ierr

		zeta_eval = zeta
		theta_eval = theta
		call spline_modulo(spline_B(isurf),theta_eval,zeta_eval)

		call ezspline_interp(spline_B(isurf),theta_eval,zeta_eval,&
			BB,ierr)
		compute_geometry_spline(B_index) = BB
		call ezspline_error(ierr)

		call ezspline_interp(spline_dBdtheta(isurf),theta_eval,zeta_eval,&
				dBBdtheta,ierr)
		compute_geometry_spline(dBdtheta_index) = dBBdtheta
		call ezspline_error(ierr)

		if (geometry_option == 2) then ! pest
			call ezspline_interp(spline_Bdotgradzeta(isurf),theta_eval,zeta_eval,&
				compute_geometry_spline(Bdotgradzeta_index),ierr)
			call ezspline_error(ierr)
			call ezspline_interp(spline_dBdotgradzetadtheta(isurf),theta_eval,zeta_eval,&
				compute_geometry_spline(dBdotgradzetadtheta_index),ierr)
			call ezspline_error(ierr)
			if (output_P_tensor) then
				call ezspline_interp(spline_d2Bdotgradzetadtheta2(isurf),theta_eval,zeta_eval,&
					compute_geometry_spline(d2Bdotgradzetadtheta2_index),ierr)
			end if
		else ! boozer
			sign_Bdotgradzeta = sign(one,Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf))
			compute_geometry_spline(Bdotgradzeta_index) = abs((BB**2) &
				/(Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf)))
			compute_geometry_spline(dBdotgradzetadtheta_index) = 2.0*BB*dBBdtheta &
				/(Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf))*sign_Bdotgradzeta
			if (output_P_tensor) then
				compute_geometry_spline(d2Bdotgradzetadtheta2_index) = 2.0*(dBBdtheta**2 + BB*d2BBdtheta2) &
					/(Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf))*sign_Bdotgradzeta
			end if
		end if

	end function compute_geometry_spline

! ===================================================
! Function compute_d2Bdtheta2_spline
!
! This function computes the local d2Bdtheta2
! using the EzSpline interface. As this is only needed
! for the bounce integrals needed to compute the P
! tensor, it has been separated from compute_geometry_spline.
!
! Input:
! isurf		index in ssurf for calculation.
! theta		theta for calculation
! zeta		zeta for calculation.
!
! Output:
! compute_d2Bdtheta2_spline		Interpolated value of
!															d2Bdtheta2
!
! ===================================================
	real(dp) function compute_d2Bdtheta2_spline(isurf,theta,zeta)

		integer, intent(in) :: isurf
		real(dp), intent(in) :: theta,zeta
		real(dp) :: theta_eval, zeta_eval
		integer :: ierr

		zeta_eval = zeta
		theta_eval = theta
		call spline_modulo(spline_B(isurf),theta_eval,zeta_eval)

		call ezspline_interp(spline_d2Bdtheta2(isurf),theta_eval,zeta_eval,&
			compute_d2Bdtheta2_spline,ierr)
		call ezspline_error(ierr)

	end function compute_d2Bdtheta2_spline


end module splines_mod
