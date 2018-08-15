module flintegrate_mod

	use stel_constants
	use stel_kinds
	use omp_lib
	use ezspline
	use splines_mod, only: spline_B, spline_dBdtheta, spline_dBdzeta
	use geometry_mod, only: iota, compute_B, nfp, Boozer_I, Boozer_G
	use constants_mod, only: integrand_length, dKdalpha_index, H_index, I_index, J_index

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
		use diagnostics_mod, only: J_invariant, dKdalpha, I_bounce_integral, one_over_nu_metric_before_integral, K_bounce_integral, nemov_metric_before_integral, H_bounce_integral, nclass
		use input_mod, only: verbose, Delta_zeta, max_search_in_zeta, nintegral, nwell, nalpha

		implicit none

		integer, intent(in) :: isurf, ilambda
		integer :: ialpha, k_well, izeta
		real(dp) :: lambda_scaled, zeta, theta
		real(dp) :: modB_left, theta_left, zeta_left, leftmost_allowed_zeta, modB
		real(dp) :: zeta_at_left_bounce_point, zeta_at_right_bounce_point
		logical :: failed
		real(dp), dimension(nintegral) :: zeta_integral
		real(dp) :: dKdalpha_integral, I_integral, J_integral, H_integral
		real(dp) :: zeta_right
		real(dp), dimension(integrand_length) :: integral
		real(dp) :: this_lambda, alpha0

		lambda_scaled = lambdas(ilambda)
		this_lambda = 1/max_B(isurf) + lambda_scaled*(1/min_B(isurf) - 1/max_B(isurf))

		!$OMP PARALLEL
		!$OMP MASTER
		if (verbose) print *,"  Number of OpenMP threads:",omp_get_num_threads()
		!$OMP END MASTER
		!$OMP DO PRIVATE(ialpha,alpha0,zeta,zeta_left,theta_left,modB_left,k_well,theta,modB,leftmost_allowed_zeta,failed,zeta_at_left_bounce_point,zeta_at_right_bounce_point,izeta,zeta_integral,dKdalpha_integral,I_integral,J_integral,H_integral,integral,zeta_right)
		do ialpha=1,nalpha
			alpha0 = alphas(ialpha)
			zeta = 0
			zeta_left = -Delta_zeta
			theta_left = alpha0 + iota(isurf)*zeta_left
			modB_left = compute_B(isurf,theta_left,zeta_left)

			k_well = 0
			do while (zeta .le. twopi/nfp)

				theta = alpha0 + iota(isurf)*zeta
				modB = compute_B(isurf,theta,zeta)

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
						modB = compute_B(isurf,theta,zeta)
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

					! Split up integral into nintegral pieces
					do izeta = 1,nintegral
						zeta_integral(izeta) = ((izeta-1.0)/(nintegral-1.0))*(zeta_at_right_bounce_point-zeta_at_left_bounce_point) + zeta_at_left_bounce_point
					end do

					dKdalpha_integral = 0
					I_integral = 0
					J_integral = 0
					H_integral = 0
					do izeta=1,(nintegral-1)
						zeta_left = zeta_integral(izeta)
						zeta_right = zeta_integral(izeta+1)
						call rk4_integrate(isurf,alpha0,this_lambda,zeta_left,zeta_right,integral)
						dKdalpha_integral = dKdalpha_integral + integral(dKdalpha_index)
						H_integral = H_integral + integral(H_index)
						I_integral = I_integral + integral(I_index)
						J_integral = J_integral + integral(J_index)
					end do
					if (k_well < nwell) then
						J_invariant(isurf,ilambda,ialpha,k_well) = J_integral
						dKdalpha(isurf,ilambda,ialpha,k_well) = dKdalpha_integral
						I_bounce_integral(isurf,ilambda,ialpha,k_well) = I_integral
						H_bounce_integral(isurf,ilambda,ialpha,k_well) = H_integral
					else
						J_invariant(isurf,ilambda,ialpha,nwell) = J_invariant(isurf,ilambda,ialpha,nwell) + J_integral
						dKdalpha(isurf,ilambda,ialpha,nwell) = dKdalpha(isurf,ilambda,ialpha,nwell) + dKdalpha_integral
						I_bounce_integral(isurf,ilambda,ialpha,nwell) = I_bounce_integral(isurf,ilambda,ialpha,nwell) + I_integral
						H_bounce_integral(isurf,ilambda,ialpha,k_well) = H_bounce_integral(isurf,ilambda,ialpha,nwell) + H_integral
					end if
				end if ! (1-lambda*modB) .ge. 0
				modB_left = modB
				zeta = zeta + Delta_zeta
			end do ! zeta
			nclass(isurf,ilambda,ialpha) = k_well
		end do ! ialpha
		!$OMP END DO
		!$OMP END PARALLEL

	end subroutine flintegrate

 ! ===================================================
 ! Subroutine bounce_integrand
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
 ! zeta			Value of zeta for evaluation of integrand
 !
 ! Outputs:
 ! result		Value of integrand
 !
 ! ===================================================
	subroutine bounce_integrand(isurf,alpha0,ilambda,zeta,result)

		integer, intent(in) :: isurf
		real(dp), intent(in) :: alpha0, ilambda, zeta
		real(dp), dimension(integrand_length) :: result

		real(dp) :: BB, dBBdtheta, dBBdzeta, theta, radicand, zeta_eval
		integer :: ierr

		theta = alpha0 + iota(isurf)*zeta
		zeta_eval = zeta

		call ezspline_modulo(spline_B(isurf),theta,zeta_eval,ierr)
		call ezspline_error(ierr)
		call ezspline_interp(spline_B(isurf),theta,zeta_eval,BB,ierr)
		call ezspline_error(ierr)
		call ezspline_interp(spline_dBdtheta(isurf),theta,zeta_eval,dBBdtheta,ierr)
		call ezspline_error(ierr)
		call ezspline_interp(spline_dBdzeta(isurf),theta,zeta_eval,dBBdzeta,ierr)
		call ezspline_error(ierr)

		radicand = max(1-ilambda*BB,0.0)

		result(dKdalpha_index) = -2*(Boozer_G(isurf)+iota(isurf)*Boozer_I(isurf))*dBBdtheta*(3*sqrt(radicand)*ilambda/(2*BB) &
			+ sqrt(radicand)**3/BB**2)
		result(H_index) = 2*(sqrt(radicand)**3 + 3*sqrt(radicand))*(Boozer_I(isurf)*dBBdzeta &
			- Boozer_G(isurf)*dBBdtheta)/(BB**3)
		result(I_index) = 2*(Boozer_G(isurf)+iota(isurf)*Boozer_I(isurf)) &
			*sqrt(radicand)/BB**2
		result(J_index) = 2*(Boozer_G(isurf)+iota(isurf)*Boozer_I(isurf)) &
			*sqrt(radicand)/BB

	end subroutine bounce_integrand

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
!
! Outputs:
! integral		Value of integral from zeta_left to zeta_right
!
! ===================================================
	subroutine rk4_integrate(this_surf,alpha0,this_lambda,zeta_left,zeta_right,integral)

		use stel_kinds
		use constants_mod

		implicit none

		integer, intent(in) :: this_surf
		real(dp), intent(in) :: alpha0, this_lambda
		real(dp), intent(in) :: zeta_left, zeta_right
		real(dp), intent(out), dimension(integrand_length) :: integral
		real(dp) :: h, hh, h6, xh, x
		real(dp), dimension(integrand_length) :: y, yt, dydx, dyt, dym

		x = zeta_left
		h = zeta_right - zeta_left
		hh = h*0.5
		h6 = h/6.0
		xh = x + hh
		y = 0

		call bounce_integrand(this_surf,alpha0,this_lambda,x,dydx)
		yt = y + hh*dydx

		call bounce_integrand(this_surf,alpha0,this_lambda,xh,dyt)
		yt = y + hh*dyt

		call bounce_integrand(this_surf,alpha0,this_lambda,xh,dym)
		yt = y + h*dym
		dym = dym + dyt

		call bounce_integrand(this_surf,alpha0,this_lambda,x+h,dyt)

		integral = y + h6*(dydx+dyt+2.0*dym)

	end subroutine rk4_integrate

end module flintegrate_mod
