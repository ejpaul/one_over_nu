module diagnostics_mod

	use stel_kinds

	! Bounce integrals
	real(dp), dimension(:,:,:,:), allocatable :: J_invariant, dKdalpha, I_bounce_integral, one_over_nu_metric_before_integral, K_bounce_integral, nemov_metric_before_integral, H_bounce_integral
	! Number of classes for each grid point
	integer, dimension(:,:,:), allocatable :: nclass
	! Output metrics
	real(dp), dimension(:), allocatable :: one_over_nu_metric, nemov_metric
	! Paritlce flux - dimension T m^{-1} s^{-1}
	real(dp), dimension(:), allocatable :: particleflux

	contains

 ! ===================================================
 ! Subroutine init_diagnostics
 !
 ! This subroutine allocates and initializes the
 ! the arrays used for diagnostics.
 !
 ! ===================================================
	subroutine init_diagnostics()

		use input_mod, only: nsurf, nlambda, nalpha, nwell

		implicit none

		allocate(one_over_nu_metric_before_integral(nsurf,&
			nlambda,nalpha,nwell))
		allocate(nemov_metric_before_integral(nsurf,nlambda,&
			nalpha,nwell))
		one_over_nu_metric_before_integral = 0
		nemov_metric_before_integral = 0

		allocate(one_over_nu_metric(nsurf))
		allocate(nemov_metric(nsurf))
		one_over_nu_metric = 0
		nemov_metric = 0

		allocate(J_invariant(nsurf,nlambda,nalpha,nwell))
		allocate(I_bounce_integral(nsurf,nlambda,nalpha,nwell))
		allocate(H_bounce_integral(nsurf,nlambda,nalpha,nwell))
		allocate(dKdalpha(nsurf,nlambda,nalpha,nwell))
		allocate(nclass(nsurf,nlambda,nalpha))
		J_invariant = 0
		I_bounce_integral = 0
		H_bounce_integral = 0
		dKdalpha = 0
		nclass = 0

	end subroutine init_diagnostics

 ! ===================================================
 ! Subroutine compute_diagnostics
 !
 ! This subroutine computes the integrals over alpha
 ! and lambda to compute output diagnostics.
 !
 ! Inputs:
 ! isurf		The index in the s_surf array corresponding
 ! 					to the surface to be calculated
 !
 ! ===================================================
	subroutine compute_diagnostics(isurf)

		use extrema_mod, only: max_B, min_B
		use grids_mod, only: dlambda, dalpha, lambdas
		use input_mod, only: nlambda, nalpha, output_particle_flux

		implicit none

		integer, intent(in) :: isurf
		real(dp) :: lambda, lambda_scaled
		integer :: ialpha, ilambda, iclass

		! Sum over particle classes
		do ilambda=1,Nlambda
			do ialpha=1,Nalpha
				do iclass=1,nclass(isurf,ilambda,ialpha)
					lambda_scaled = lambdas(ilambda)
					lambda = 1.0/max_B(isurf) + &
				(1.0/min_B(isurf) - 1.0/max_B(isurf)) * lambda_scaled
					one_over_nu_metric_before_integral(isurf,ilambda,ialpha,iclass) = &
						(dKdalpha(isurf,ilambda,ialpha,iclass)**2) &
						/I_bounce_integral(isurf,ilambda,ialpha,iclass) &
						*(1.0/lambda)
					nemov_metric_before_integral(isurf,ilambda,ialpha,iclass) = &
						(H_bounce_integral(isurf,ilambda,ialpha,iclass)**2) &
						/I_bounce_integral(isurf,ilambda,ialpha,iclass) &
						*(1.0/lambda)
				end do
			end do
		end do

		! Trapezoid rule for lambda integral
		one_over_nu_metric_before_integral(isurf,1,:,:) = &
			0.5*one_over_nu_metric_before_integral(isurf,1,:,:)
		one_over_nu_metric_before_integral(isurf,nlambda,:,:) = &
			0.5*one_over_nu_metric_before_integral(isurf,nlambda,:,:)
		nemov_metric_before_integral(isurf,1,:,:) = 0.5*nemov_metric_before_integral(isurf,1,:,:)
		nemov_metric_before_integral(isurf,nlambda,:,:) = &
			0.5*nemov_metric_before_integral(isurf,nlambda,:,:)

		! Compute integrals
		one_over_nu_metric(isurf) = &
			sum(one_over_nu_metric_before_integral(isurf,:,:,:))*dalpha*dlambda(isurf)
		nemov_metric(isurf) = sum(nemov_metric_before_integral(isurf,:,:,:))*dalpha*dlambda(isurf)

		if (output_particle_flux > 0) then
			call compute_particle_flux()
		end if

	end subroutine compute_diagnostics

! ===================================================
! Subroutine compute_particle flux
!
! This subroutine computes the particle flux
! using user-specified density and temperature
! profiles. This is only called if
! output_particle_flux_option = 1.
!
! ===================================================
	subroutine compute_particle_flux

		use input_mod
		use splines_mod
		use geometry_mod
		use stel_constants
		use constants_mod

		implicit none

		real(dp) :: lnlambda, nuhat, v_t1, vprime, prefactor, H1, H2
		integer :: isurf

		allocate(particleFlux(nsurf))

		do isurf = 1, nsurf
			if (q_e(1)<0 .and. q_e(2)<0) then
				! Defined in NRL formulary - electron-electron self collisions
				lnlambda = 23.5 - 0.5*log(n_m3(isurf)*1.e-6) + 1.25*log(T_ev(isurf)) &
					- sqrt(1.0e-5 + 0.0625*(log(T_ev(isurf))-2)**2)
			else if (q_e(1)<0 .and. q_e(2)>0) then
				! Defined in NRL formulary - electron-ion collisions
				lnlambda = 24 - 0.5*log(n_m3(isurf)*(1.0e-6)) + log(T_ev(isurf))
			else
				! Defined in NRL formular - ion-ion self collisions
				lnlambda = 23 - log(q_e(1)*q_e(2)/T_ev(isurf)) &
					- 0.5*log(2.0*n_m3(isurf)*(1.0e-6)*(q_e(1)**2)/T_ev(isurf))
			end if

			! Compute thermal speed for primary species
			v_t1 = sqrt(2*T_ev(isurf)*e/(m_kg))
			! Compute collision frequency
			nuhat = n_m3(isurf)*q_e(2)*(q_e(1)**3)*(e**4)*lnlambda/ &
				(4*pi*(epsilon0**2)*(m_kg**2)*(v_t1**3))

			! Compute vprime using spline grid
			vprime = (Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf)) &
				*sum(1.0/B_for_spline(isurf,:,:)**2)*dtheta_spline*dzeta_spline

			prefactor = -4*n_m3(isurf)*(T_eV(isurf)**2)/(9*vprime*nuhat*sqrt(pi))
			H1 = prefactor*6.770
			H2 = prefactor*29.45

			particleFlux(isurf) = one_over_nu_metric(isurf) &
				*(H1*dlnndpsi(isurf) + (H2 - 1.5*H1)*dlnTdpsi(isurf))
		end do

	end subroutine compute_particle_flux

end module diagnostics_mod
