module extrema_mod

	use stel_kinds
	use input_mod, only: niter_newton, tol_newton

	! Input parameters
	integer :: ntheta = 100 ! Number of points in equally-spaced theta grid for finding (rough) extrema
	integer :: nzeta = 100 ! Number of points in equally-spaced theta grid for finding (rough) extrema
	! Grids
	real(dp), dimension(:), allocatable :: thetas, zetas
	real(dp), dimension(:,:,:), allocatable :: B
	! Extrema of modB
	real(dp), dimension(:), allocatable :: min_B, max_B

	contains

 ! ===================================================
 ! Subroutine find_extrema
 !
 ! This subroutine finds the min and max B, first
 ! approximately on the ntheta-nzeta grid, then
 ! exactly using a Newton solve.
 !
 ! ===================================================
	subroutine find_extrema()

		use grids_mod, only: dlambda
		use geometry_mod, only: compute_B, nfp
		use input_mod, only: verbose, nsurf
		use grids_mod, only: dlambda
		use stel_constants

		implicit none

		integer, dimension(2) :: result
		real(dp) :: theta_min, zeta_min, theta_max, zeta_max
		integer :: theta_min_coarse, zeta_min_coarse, theta_max_coarse, zeta_max_coarse, isurf
		integer :: itheta, izeta

		allocate(min_B(nsurf))
		allocate(max_B(nsurf))
		allocate(thetas(ntheta))
		allocate(zetas(nzeta))
		allocate(B(nsurf,ntheta,nzeta))

		! Construct theta-zeta grid for finding approximate extrema
		do izeta=1,nzeta
			zetas(izeta) = (twopi/nfp)*(izeta-1.0)/(nzeta)
		end do

		do itheta=1,ntheta
			thetas(itheta) = twopi*(itheta-1.0)/(ntheta)
		end do

		! Compute B on theta-zeta grid
		B = 0
		do isurf=1,nsurf
			do itheta=1,ntheta
				do izeta=1,nzeta
					B(isurf,itheta,izeta) = compute_B(isurf,thetas(itheta),zetas(izeta))
				end do
			end do
		end do

		do isurf=1,nsurf

			! Find approximate minima on theta-zeta grid
			result = minloc(B(isurf,:,:))
			theta_min_coarse = result(1)
			zeta_min_coarse = result(2)
			if (verbose) then
				print *,"B_min_coarse: ", compute_B(isurf,thetas(theta_min_coarse),zetas(zeta_min_coarse))
			end if

			! Find approximate maxima on theta-zeta grid
			result = maxloc(B(isurf,:,:))
			theta_max_coarse = result(1)
			zeta_max_coarse = result(2)
			if (verbose) then
				print *,"B_max_coarse: ", compute_B(isurf,thetas(theta_max_coarse),zetas(zeta_max_coarse))
			end if

			! Find exact extrema using newton solve
			call newton_opt(isurf,thetas(theta_min_coarse),zetas(zeta_min_coarse),&
				theta_min,zeta_min,thetas(theta_max_coarse),zetas(zeta_max_coarse),theta_max,zeta_max)

			min_B(isurf) = compute_B(isurf,theta_min,zeta_min)
			if (verbose) then
				print *,"B_min: ", min_B(isurf)
			end if

			max_B(isurf) = compute_B(isurf,theta_max,zeta_max)
			if (verbose) then
				print *,"B_max: ", max_B(isurf)
			end if

			dlambda(isurf) = dlambda(isurf)*(1.0/max_B(isurf) - 1.0/min_B(isurf))

		end do ! isurf

	end subroutine find_extrema

! ===================================================
! Subroutine newton_opt
!
! This subroutine finds the min and max of B
! using a Newton method.
!
! Inputs:
! isurf						index of s_surf for calculation
! theta_min_rough	Theta corresponding to the approximate
!									minima
! zeta_min_rough  Zeta corresponding to the approximate
!									minima
!	theta_max_rough Theta corresponding to the approximate
!									maxima
! zeta_max_rough  Zeta corresponding to the approximate
!									maxima
!
!	Outputs:
!	theta_min				Theta corresponding to minima (exact)
! zeta_min				Zeta corresponding to minima (exact)
!	theta_max				Theta corresponding to maxima (exact)
! zeta_max				Zeta corresponding to maxima (exact)
!
! ===================================================
	subroutine newton_opt(isurf,theta_min_rough,zeta_min_rough,theta_min,&
			zeta_min,theta_max_rough,zeta_max_rough,theta_max,zeta_max)

		use geometry_mod, only: compute_geometry, xm, xn, bmnc
		use constants_mod

		implicit none

		real(dp), intent(in) :: theta_min_rough, zeta_min_rough, theta_max_rough, zeta_max_rough
		real(dp), intent(out) :: theta_min, zeta_min, theta_max, zeta_max
		integer, intent(in) :: isurf

		real(dp) :: dBdtheta,dBdzeta,dB2dtheta2,dB2dzeta2,dB2dthetadzeta
		real(dp) :: minus_dBdtheta,minus_dBdzeta,minus_dB2dtheta2,minus_dB2dzeta2,minus_dB2dthetadzeta

		real(dp) :: theta_old, zeta_old, det, theta_new, zeta_new
		real(dp) :: geometry(geometry_length)
		integer :: iter, code

		theta_old = theta_min_rough
		zeta_old = zeta_min_rough
		code = -1
		do iter=1,Niter_Newton
			! Compute gradient and Hessian matrix
			geometry = compute_geometry(isurf,theta_old,zeta_old)
			dBdtheta = geometry(dBdtheta_index)
			dBdzeta = geometry(dBdzeta_index)
			dB2dtheta2 = -dot_product(bmnc(isurf,:),(xm**2)*cos(xm*theta_old-xn*zeta_old))
			dB2dzeta2 = -dot_product(bmnc(isurf,:),(xn**2)*cos(xm*theta_old-xn*zeta_old))
			dB2dthetadzeta = dot_product(bmnc(isurf,:),xm*xn*cos(xm*theta_old-xn*zeta_old))
			det = abs(dB2dtheta2*dB2dzeta2-dB2dthetadzeta**2)
			! Take Newton step
			theta_new = theta_old - (dB2dzeta2*dBdtheta - dB2dthetadzeta*dBdzeta)/det
			zeta_new = zeta_old - (dB2dtheta2*dBdzeta - dB2dthetadzeta*dBdtheta)/det
			if (abs(theta_new - theta_old) < tol_newton .and. abs(zeta_new - zeta_old) < tol_newton) then
				code = 0
				exit
			end if
			zeta_old = zeta_new
			theta_old = theta_new
		end do
		if (code < 0) then
			print *, 'Newton iteration for min_B did not converge!'
			stop
		end if
		zeta_min = zeta_new
		theta_min = theta_new

		theta_old = theta_max_rough
		zeta_old = zeta_max_rough
		code = -1
		do iter=1,Niter_Newton
			! Compute gradient and Hessian matrix
			geometry = compute_geometry(isurf,theta_old,zeta_old)
			minus_dBdtheta = -geometry(dBdtheta_index)
			minus_dBdzeta = -geometry(dBdzeta_index)
			minus_dB2dtheta2 = dot_product(bmnc(isurf,:),(xm**2)*cos(xm*theta_old-xn*zeta_old))
			minus_dB2dzeta2 = dot_product(bmnc(isurf,:),(xn**2)*cos(xm*theta_old-xn*zeta_old))
			minus_dB2dthetadzeta = -dot_product(bmnc(isurf,:),xm*xn*cos(xm*theta_old-xn*zeta_old))
			det = abs(minus_dB2dtheta2*minus_dB2dzeta2-minus_dB2dthetadzeta**2)
			! Take Newton step
			theta_new = theta_old - (minus_dB2dzeta2*minus_dBdtheta - minus_dB2dthetadzeta*minus_dBdzeta)/det
			zeta_new = zeta_old - (minus_dB2dtheta2*minus_dBdzeta - minus_dB2dthetadzeta*minus_dBdtheta)/det
			if (abs(theta_new - theta_old) < tol_newton .and. abs(zeta_new - zeta_old) < tol_newton) then
				code = 0
				exit
			end if
			zeta_old = zeta_new
			theta_old = theta_new
		end do
		if (code < 0) then
			print *, 'Newton iteration for max_B did not converge!'
			stop
		end if
		zeta_max = zeta_new
		theta_max = theta_new

	end subroutine newton_opt

end module extrema_mod
