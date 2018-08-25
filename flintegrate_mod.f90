module flintegrate_mod

	use stel_constants
	use stel_kinds
	use omp_lib
	use ezspline
	use splines_mod, only: compute_geometry_spline, compute_B_spline
	use geometry_mod, only: iota, nfp
	use constants_mod

	implicit none

	contains

! ===================================================
! Subroutine flintegrate
!
! This subroutine performs moves a long a field line
! to find the bounce points and perform the required
! bounce integrals.
!
! Inputs:
!	isurf		integer in ssurf corresponding to surface
!					for calculation
! ilambda	integer in lambdas corresponding to
!					(scaled) lambda for calculation
!
! ===================================================
	subroutine flintegrate(isurf,ilambda)

		use grids_mod, only: alphas, lambdas
		use extrema_mod, only: min_B, max_B
		use diagnostics_mod, only: J_invariant, dKdalpha, I_bounce_integral, &
			one_over_nu_metric_before_integral, K_bounce_integral,nclass, d2Kdalpha2, dIdalpha
		use input_mod, only: verbose, Delta_zeta, max_search_in_zeta, nintegral, nwell, nalpha, output_P_tensor, output_J

		implicit none

		integer, intent(in) :: isurf, ilambda
		integer :: ialpha, k_well, izeta
		real(dp) :: lambda_scaled, zeta, theta
		real(dp) :: modB_left, theta_left, zeta_left, leftmost_allowed_zeta, modB
		real(dp) :: zeta_at_left_bounce_point, zeta_at_right_bounce_point
		logical :: failed
		real(dp), dimension(integrand_length) :: integral
		real(dp) :: alpha0, this_lambda

		lambda_scaled = lambdas(ilambda)
		this_lambda = 1/max_B(isurf) + lambda_scaled*(1/min_B(isurf) - 1/max_B(isurf))

		!$OMP PARALLEL
		!$OMP MASTER
		if (verbose) print *,"  Number of OpenMP threads:",omp_get_num_threads()
		!$OMP END MASTER
		!$OMP DO PRIVATE(ialpha,alpha0,zeta,zeta_left,theta_left,modB_left,k_well,theta,modB,leftmost_allowed_zeta,failed,zeta_at_left_bounce_point,zeta_at_right_bounce_point,izeta,integral)
		do ialpha=1,nalpha
			alpha0 = alphas(ialpha)
			zeta = 0
			zeta_left = -Delta_zeta
			theta_left = alpha0 + iota(isurf)*zeta_left
			modB_left = compute_B_spline(isurf,theta_left,zeta_left)

			k_well = 0
			do while (zeta .le. twopi/nfp)

				theta = alpha0 + iota(isurf)*zeta
				modB = compute_B_spline(isurf,theta,zeta)

				if ((1-this_lambda*modB) .ge. 0) then
					leftmost_allowed_zeta = zeta
					if (((1-this_lambda*modB_left)>0) .or. ((leftmost_allowed_zeta - Delta_zeta)<0)) then
						modB_left = modB
						zeta = zeta + Delta_zeta
						cycle
					end if

					! Accurately solve for left bounce point
					call find_bounce_point(isurf,leftmost_allowed_zeta-Delta_zeta,leftmost_allowed_zeta,alpha0,this_lambda,zeta_at_left_bounce_point)

					! Step to right to roughly locate right bounce point
					failed = .true.

					do while (zeta .le. max_search_in_zeta)
						! step to the right
						zeta = zeta + Delta_zeta
						! Check whether this point is forbidden
						theta = alpha0 + iota(isurf)* zeta
						modB = compute_B_spline(isurf,theta,zeta)

						if ((1 - this_lambda*modB) < 0) then
							failed = .false.
							exit
						end if
					end do
					if (failed) then
						print *,"Could not find right bounce point!"
						print *,"ialpha: ", ialpha
						print *,"ilambda: ", ilambda
						cycle
					end if

					! Accurately solve for right bounce point
					call find_bounce_point(isurf,zeta-Delta_zeta,zeta,alpha0,this_lambda,zeta_at_right_bounce_point)

					k_well = k_well + 1

					! P tensor requires computing bounce integrals with sqrt(1-lambda*B) in denominator
					! so integration performed with variable transformation.
					if (output_P_tensor) then
						call bounce_integrate_remapped(isurf,alpha0,this_lambda,zeta_at_left_bounce_point,&
							zeta_at_right_bounce_point,integral)
					else
						call bounce_integrate(isurf,alpha0,this_lambda,zeta_at_left_bounce_point,&
							zeta_at_right_bounce_point,integral)
					end if

					if (k_well < nwell) then
						if (output_J) then
							J_invariant(isurf,ilambda,ialpha,k_well) = integral(J_index)
						end if
						dKdalpha(isurf,ilambda,ialpha,k_well) = integral(dKdalpha_index)
						I_bounce_integral(isurf,ilambda,ialpha,k_well) = integral(I_index)
						if (output_P_tensor) then
							d2Kdalpha2(isurf,ilambda,ialpha,k_well) = integral(d2Kdalpha2_index)
							dIdalpha(isurf,ilambda,ialpha,k_well) = integral(dIdalpha_index)
						end if
					else
						if (output_J) then
							J_invariant(isurf,ilambda,ialpha,nwell) = J_invariant(isurf,ilambda,ialpha,nwell) + integral(J_index)
						end if
						dKdalpha(isurf,ilambda,ialpha,nwell) = dKdalpha(isurf,ilambda,ialpha,nwell) + integral(dKdalpha_index)
						I_bounce_integral(isurf,ilambda,ialpha,nwell) = I_bounce_integral(isurf,ilambda,ialpha,nwell) + integral(I_index)
						if (output_P_tensor) then
							d2Kdalpha2(isurf,ilambda,ialpha,nwell) = d2Kdalpha2(isurf,ilambda,ialpha,nwell) + integral(d2Kdalpha2_index)
							dIdalpha(isurf,ilambda,ialpha,nwell) = dIdalpha(isurf,ilambda,ialpha,nwell) + integral(dIdalpha_index)
						end if
					end if
				end if
				modB_left = modB
				zeta = zeta + Delta_zeta
			end do ! zeta
			nclass(isurf,ilambda,ialpha) = k_well
		end do ! ialpha
		!$OMP END DO
		!$OMP END PARALLEL

	end subroutine flintegrate

! ===================================================
! Subroutine bounce_integrate_remapped
!
! This subroutine is used to compute bounce integrals
! that contain an integrable singularity. The domain
! is split into two, where the variable t = sqrt(zeta - zeta_b)
! is used on the left side and t = sqrt(zeta_b - zeta)
! on the right side. Nintegral rk4 steps are performed
! on each side to compute the bounce integral.
!
! Input:
! isurf		integer in ssurf corresponding to surface
!					for calculation
! alpha0
! ilambda  integer in lambdas corresponding to
!					(scaled) lambda for calculation
! zeta_at_left_bounce_point		value of zeta for left bounce point
! zeta_at_right_bounce_point  value of zeta for right bounce point
!
! Output:
!	integral	Value of bounce integral. This is an array
! 					of length integrand_length where the array
!						can be indexed using parameters defined in
!						constants_mod.f90
!
! ===================================================
	subroutine bounce_integrate_remapped(isurf,alpha0,this_lambda,zeta_at_left_bounce_point,&
		zeta_at_right_bounce_point,integral)

		use input_mod, only: nintegral

		integer, intent(in) :: isurf
		real(dp), intent(in) :: alpha0, this_lambda, zeta_at_left_bounce_point, zeta_at_right_bounce_point
		real(dp), intent(out), dimension(integrand_length) :: integral
		real(dp) :: midpoint, zeta_b, t_left, t_right
		real(dp), dimension(nintegral+1) :: t_integral
		real(dp), dimension(integrand_length) :: this_integral, integral_left, integral_right
		integer :: it, which_integral

		midpoint = zeta_at_left_bounce_point + 0.5*(zeta_at_right_bounce_point - zeta_at_left_bounce_point)

		! Perform integration from zeta_at_left_bounce_point to midpoint
		zeta_b = zeta_at_left_bounce_point
		which_integral = -1
		t_left = 0
		t_right = sqrt(midpoint-zeta_b)
		! Split up integral into Nintegral pieces
		do it=1,Nintegral+1
			t_integral(it) = ((it-1.0)/(nintegral))*(t_right-t_left) + t_left
		end do

		integral_left = 0.0
		do it = 1,nintegral
			t_left = t_integral(it)
			t_right = t_integral(it+1)
			call rk4_integrate(isurf,alpha0,this_lambda,t_left,t_right,this_integral,&
				zeta_b,which_integral)
			integral_left = integral_left + this_integral
		end do

		! Perform integration from midpoint to zeta_at_right_bounce_point
		zeta_b = zeta_at_right_bounce_point
		which_integral = 1
		t_left = sqrt(zeta_b-midpoint)
		t_right = 0
		! Split up integral into Nintegral pieces
		do it=1,Nintegral+1
			t_integral(it) = ((it-1.0)/(nintegral))*(t_right-t_left) + t_left
		end do
		integral_right = 0.0
		do it = 1,nintegral
			t_left = t_integral(it)
			t_right = t_integral(it+1)
			call rk4_integrate(isurf,alpha0,this_lambda,t_left,t_right,this_integral,&
				zeta_b,which_integral)
			integral_right = integral_right + this_integral
		end do

		integral = integral_left + integral_right

	end subroutine bounce_integrate_remapped

! ===================================================
! Subroutine bounce_integrate
!
! This subroutine is used to compute bounce integrals
! that don't contain an integrable singularity. Nintegral
! rk4 steps are performed to compute the bounce integral.
!
! Input:
! isurf		integer in ssurf corresponding to surface
!					for calculation
! alpha0
! ilambda  integer in lambdas corresponding to
!					(scaled) lambda for calculation
! zeta_at_left_bounce_point		value of zeta for left bounce point
! zeta_at_right_bounce_point  value of zeta for right bounce point
!
! Output:
!	integral	Value of bounce integral. This is an array
! 					of length integrand_length where the array
!						can be indexed using parameters defined in
!						constants_mod.f90
!
! ===================================================
	subroutine bounce_integrate(isurf,alpha0,this_lambda,zeta_at_left_bounce_point,&
		zeta_at_right_bounce_point,integral)

			use input_mod, only: nintegral

			implicit none

			integer, intent(in) :: isurf
			real(dp), intent(in) :: alpha0, this_lambda, zeta_at_left_bounce_point,zeta_at_right_bounce_point
			real(dp), dimension(integrand_length), intent(out) :: integral
			real(dp), dimension(integrand_length) :: this_integral
			integer :: izeta
			real(dp) :: zeta_left, zeta_right
			real(dp), dimension(nintegral+1) :: zeta_integral
			real(dp) :: zeta_b = -1
			integer :: which_integral = 0

			! Split up integral into nintegral pieces
			do izeta = 1,nintegral+1
				zeta_integral(izeta) = ((izeta-1.0)/(nintegral))*(zeta_at_right_bounce_point-zeta_at_left_bounce_point) + zeta_at_left_bounce_point
			end do
			integral = 0.0
			do izeta = 1,(nintegral)
				zeta_left = zeta_integral(izeta)
				zeta_right = zeta_integral(izeta+1)
				call rk4_integrate(isurf,alpha0,this_lambda,zeta_left,zeta_right,&
					this_integral,zeta_b,which_integral)
				integral = integral + this_integral
			end do

	end subroutine bounce_integrate

 ! ===================================================
 ! Function bounce_integrand
 !
 ! This subroutine computes the integrands for
 ! bounce integration. This is called by the
 ! Runge-Kutta integration routine.
 !
 ! Inputs:
 ! isurf		integer in ssurf corresponding to surface
 !					for calculation
 ! alpha0
 ! ilambda  integer in lambdas corresponding to
 !					(scaled) lambda for calculation
 ! t			  If which_integrand = 0, Value of zeta for
 !					evaluation of
 !					integrand. Else, value of t for evaluation of
 !					integrand.
 ! zeta_b		If which_integrand = 0, then this is the
 !					integrand of
 ! 					the zeta integral. Else, zeta_b is the bounce
 ! 					point corresponding to the remapped integral
 !					with variable t = sqrt(zeta - zeta_b) or
 !					sqrt(zeta_b - zeta)
 ! which_integrand		If which_integrand = 0, integral is
 !										performed with respect to zeta. If
 !										which_integrand = -1, integral performed
 !										with respect to t = sqrt(zeta - zeta_b). If
 !										which_integrand = 1, integral performed
 !										with respect to t = sqrt(zeta_b - zeta).
 ! Outputs:
 ! bounce_integrand		Value of integrand
 !
 ! ===================================================
	function bounce_integrand(isurf,alpha0,lambda,t,zeta_b,which_integral)

		use constants_mod
		use input_mod
		use splines_mod
		use geometry_mod, only: Boozer_G, Boozer_I, iota

		integer, intent(in) :: isurf, which_integral
		real(dp), intent(in) :: alpha0, lambda, t, zeta_b
		real(dp), dimension(integrand_length) :: bounce_integrand
		real(dp), dimension(geometry_length) :: geometry
		real(dp) :: BB, dBBdtheta, dBBdzeta, theta, radicand, d2BBdtheta2, zeta, dzetadt, &
			inv_radicand, BBdotgradzeta, dBBdotgradzetadtheta, d2BBdotgradzetadtheta2
		integer :: ierr

		if (which_integral == 0) then
			zeta = t
		! Integral from zeta_at_left_bounce_point to midpoint
		else if (which_integral == -1) then
			zeta = (t**2) + zeta_b
			dzetadt = 2.0*t
		! Integral from midpoint to zeta_right_bounce_point
		else if (which_integral == 1) then
			zeta = zeta_b - (t**2)
			dzetadt = -2.0*t
		else
			stop "Incorrect argument which_integral to bounce_integrand."
		end if

		theta = alpha0 + iota(isurf)*zeta

		geometry = compute_geometry_spline(isurf,theta,zeta)
		BB = geometry(B_index)
		dBBdtheta = geometry(dBdtheta_index)
		BBdotgradzeta = geometry(Bdotgradzeta_index)
		dBBdotgradzetadtheta = geometry(dBdotgradzetadtheta_index)

		radicand = max((1-lambda*BB),0.0)

		bounce_integrand(dKdalpha_index) = -3.0*sqrt(radicand)*lambda*dBBdtheta/BBdotgradzeta &
			 - 2*(sqrt(radicand)**3)*dBBdotgradzetadtheta/(BBdotgradzeta**2)
		bounce_integrand(I_index) = 2.0*sqrt(radicand)/BBdotgradzeta
		if (output_J) then
			bounce_integrand(J_index) = 2.0*BB*sqrt(radicand)/BBdotgradzeta
		end if
		! Compute bounce integrals for P tensor
		if (which_integral .ne. 0) then
			inv_radicand = max(1.0/(1-lambda*BB),0.0)
			d2BBdtheta2 = compute_d2Bdtheta2_spline(isurf,theta,zeta)
			d2BBdotgradzetadtheta2 = geometry(d2Bdotgradzetadtheta2_index)
			bounce_integrand(d2Kdalpha2_index) = -3.0*sqrt(radicand)*lambda*d2BBdtheta2/BBdotgradzeta &
				+ 1.5*lambda*lambda*(dBBdtheta**2)*sqrt(inv_radicand)/BBdotgradzeta &
				+ 3.0*sqrt(radicand)*lambda*dBBdtheta*dBBdotgradzetadtheta/(BBdotgradzeta**2) &
				- 2.0*(sqrt(radicand)**3)*d2BBdotgradzetadtheta2/(BBdotgradzeta**2) &
				+ 3.0*sqrt(radicand)*lambda*dBBdtheta*dBBdotgradzetadtheta/BBdotgradzeta &
				+ 4.0*(sqrt(radicand)**3)*(dBBdotgradzetadtheta**2)/(BBdotgradzeta**3)
			bounce_integrand(dIdalpha_index) = -3.0*lambda*dBBdtheta*sqrt(inv_radicand)/BBdotgradzeta &
				- 2.0*sqrt(radicand)*dBBdotgradzetadtheta/(BBdotgradzeta**2)
			bounce_integrand = bounce_integrand*dzetadt
		end if

	end function bounce_integrand

! ===================================================
! Subroutine rk4_integrate
!
! This subroutine performs integration along a field
! line from zeta_left and zeta_right using a 4th
! order Runge-Kutta method.
!
! Inputs:
! this_surf		Integer corresponding to magnetic surface for integral.
! alpha0			Field line label for integral.
! this_lambda Value of lambda for integral.
! zeta_left		Left bound of integral over zeta.
! zeta_right  Right bounc of integral over zeta.
! which_integral	Denotes if integral performed wrt zeta or t as
!									defined in bounce_integrand subroutine.
!
! Outputs:
! integral		Value of integral from zeta_left to zeta_right
!
! ===================================================
	subroutine rk4_integrate(this_surf,alpha0,this_lambda,zeta_left,&
			zeta_right,integral,zeta_b, which_integral)

		use stel_kinds
		use constants_mod

		implicit none

		integer, intent(in) :: this_surf, which_integral
		real(dp), intent(in) :: alpha0, this_lambda
		real(dp), intent(in) :: zeta_left, zeta_right, zeta_b
		real(dp), intent(out), dimension(integrand_length) :: integral
		real(dp) :: h, hh, h6, xh, x
		real(dp), dimension(integrand_length) :: y, yt, dydx, dyt, dym

		x = zeta_left
		h = zeta_right - zeta_left
		hh = h*0.5
		h6 = h/6.0
		xh = x + hh
		y = 0

		dydx = bounce_integrand(this_surf,alpha0,this_lambda,x,zeta_b,which_integral)
		yt = y + hh*dydx

		dyt = bounce_integrand(this_surf,alpha0,this_lambda,xh,zeta_b,which_integral)
		yt = y + hh*dyt

		dym = bounce_integrand(this_surf,alpha0,this_lambda,xh,zeta_b,which_integral)
		yt = y + h*dym
		dym = dym + dyt

		dyt = bounce_integrand(this_surf,alpha0,this_lambda,x+h,zeta_b,which_integral)

		integral = y + h6*(dydx+dyt+2.0*dym)

	end subroutine rk4_integrate

end module flintegrate_mod
