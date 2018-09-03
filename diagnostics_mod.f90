module diagnostics_mod

	use stel_kinds

	! Bounce integrals
	real(dp), dimension(:,:,:,:), allocatable :: J_invariant, dKdalpha, I_bounce_integral, one_over_nu_metric_before_integral, K_bounce_integral, dIdalpha, d2Kdalpha2
	! Number of classes for each grid point
	integer, dimension(:,:,:), allocatable :: nclass
	! Output metrics
	real(dp), dimension(:), allocatable :: one_over_nu_metric
	! Paritlce flux - dimension T m^{-1} s^{-1}
	real(dp), dimension(:), allocatable :: particleflux
	real(dp), dimension(:), allocatable :: energy_integral
	! P tensor for adjoint perturbed equilibrium
	real(dp), dimension(:,:,:), allocatable :: P_tensor_bb, P_tensor_I
	real(dp), dimension(:,:,:,:,:), allocatable :: P_tensor_bb_before_integral
	real(dp), dimension(:,:,:,:,:), allocatable :: P_tensor_I_before_integral

	contains

 ! ===================================================
 ! Subroutine init_diagnostics
 !
 ! This subroutine allocates and initializes the
 ! the arrays used for diagnostics.
 !
 ! ===================================================
	subroutine init_diagnostics()

		use input_mod, only: nsurf, nlambda, nalpha, nwell, output_J, &
			output_P_tensor, output_particle_flux, nzeta_spline

		implicit none

		allocate(one_over_nu_metric_before_integral(nsurf,&
			nlambda,nalpha,nwell))
		one_over_nu_metric_before_integral = 0

		allocate(one_over_nu_metric(nsurf))
		one_over_nu_metric = 0

		if (output_J) then
			allocate(J_invariant(nsurf,nlambda,nalpha,nwell))
			J_invariant = 0
		end if
		allocate(I_bounce_integral(nsurf,nlambda,nalpha,nwell))
		allocate(dKdalpha(nsurf,nlambda,nalpha,nwell))
		if (output_P_tensor) then
			allocate(d2Kdalpha2(nsurf,nlambda,nalpha,nwell))
			allocate(dIdalpha(nsurf,nlambda,nalpha,nwell))
			d2Kdalpha2 = 0
			dIdalpha = 0
		end if
		allocate(nclass(nsurf,nlambda,nalpha))

		if (output_P_tensor .or. output_particle_flux) then
			allocate(energy_integral(nsurf))
		end if

		if (output_P_tensor) then
			allocate(P_tensor_bb(nsurf,nalpha,nzeta_spline))
			allocate(P_tensor_I(nsurf,nalpha,nzeta_spline))
			allocate(P_tensor_bb_before_integral(nsurf,nalpha,nzeta_spline,nlambda,nwell))
			allocate(P_tensor_I_before_integral(nsurf,nalpha,nzeta_spline,nlambda,nwell))
			P_tensor_bb_before_integral = 0.0
			P_tensor_I_before_integral = 0.0
			P_tensor_bb_before_integral = 0.0
			P_tensor_I_before_integral = 0.0
		end if

		if (output_particle_flux) then
			allocate(particleFlux(nsurf))
		end if

		I_bounce_integral = 0
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
		use input_mod, only: nlambda, nalpha, output_particle_flux, output_p_tensor

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
				end do
			end do
		end do

		! Trapezoid rule for lambda integral
		one_over_nu_metric_before_integral(isurf,1,:,:) = &
			0.5*one_over_nu_metric_before_integral(isurf,1,:,:)
		one_over_nu_metric_before_integral(isurf,nlambda,:,:) = &
			0.5*one_over_nu_metric_before_integral(isurf,nlambda,:,:)

		! Compute integrals
		one_over_nu_metric(isurf) = &
			sum(one_over_nu_metric_before_integral(isurf,:,:,:))*dalpha*dlambda(isurf)

		if (output_particle_flux .or. output_p_tensor) then
			call compute_energy_integral(isurf)
		end if
		if (output_particle_flux) then
			call compute_particle_flux(isurf)
		end if
		if (output_p_tensor) then
			call compute_p_tensor(isurf)
		end if

	end subroutine compute_diagnostics

! ===================================================
! Subroutine compute_particle_flux
!
! This subroutine computes the particle flux
! using user-specified density and temperature
! profiles. This is only called if
! output_particle_flux_option = .true..
!
! ===================================================
	subroutine compute_particle_flux(isurf)

		use input_mod
		use splines_mod
		use geometry_mod
		use stel_constants
		use constants_mod

		implicit none

		integer, intent(in) :: isurf
		real(dp) :: lnlambda, nuhat, v_t1, vprime, prefactor, H1, H2

		! Compute vprime using spline grid
		if (geometry_option==1) then
			vprime = (Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf)) &
				*sum(1.0/B_for_spline(isurf,:,:)**2)*dtheta_spline*dzeta_spline
		else
			vprime = sum(1.0/Bdotgradzeta_for_spline(isurf,:,:))*dtheta_spline*dzeta_spline
		end if

		prefactor = -pi*(m_kg**2)/(9*vprime*(e**2)*(q_e(1)**2))

		particleFlux(isurf) = prefactor*one_over_nu_metric(isurf)*energy_integral(isurf)

	end subroutine compute_particle_flux

! ===================================================
! Subroutine compute_energy_integral
!
! This subroutine compute the energy integral needed
! for the computation of the particle flux and
! P tensor.
!
! ===================================================
	subroutine compute_energy_integral(isurf)

		use input_mod
		use stel_constants
		use constants_mod

		implicit none

		integer, intent(in) :: isurf
		real(dp) :: lnlambda, nuhat, v_t1

		if (collision_species_option == 2 .and. q_e(1)<0) then
			! Defined in NRL formulary - electron-electron self collisions
			lnlambda = 23.5 - 0.5*log(n_m3(isurf)*1.e-6) + 1.25*log(T_ev(isurf)) &
				- sqrt(1.0e-5 + 0.0625*(log(T_ev(isurf))-2)**2)
		else if (collision_species_option == 1) then
			! Defined in NRL formulary - electron-ion collisions
			lnlambda = 24 - 0.5*log(n_m3(isurf)*(1.0e-6)) + log(T_ev(isurf))
		else if (collision_species_option == 2 .and. q_e(1)>0) then
			! Defined in NRL formulary - ion-ion self collisions
			lnlambda = 23 - log(q_e(1)**2/T_ev(isurf)) &
				- 0.5*log(2.0*n_m3(isurf)*(1.0e-6)*(q_e(1)**2)/T_ev(isurf))
		else
			stop "Incorrect option in compute_particle_flux!"
		end if

		! Compute thermal speed for primary species
		v_t1 = sqrt(2*T_ev(isurf)*e/(m_kg))
		! Compute collision frequency
		if (collision_species_option==1) then
			! electron-ion collisions
			nuhat = n_m3(isurf)*q_e(2)*(-q_e(1)**3)*(e**(4))*lnlambda/ &
				(4*pi*(epsilon0**2)*(m_kg**2)*(v_t1**3))
		else
			! self collisions
			nuhat = n_m3(isurf)*(q_e(1)**4)*(e**(4))*lnlambda/ &
				(4*pi*(epsilon0**2)*(m_kg**2)*(v_t1**3))
		end if

		if (collision_species_option == 1) then
			! electron-ion collisions
			energy_integral(isurf) = 12*n_m3(isurf)*(v_t1**4)/(nuhat*(pi**(1.5))) &
				*(dlnndpsi(isurf) + 3.5*dlnTdpsi(isurf))
		else
			! self collisions
			energy_integral(isurf) = 13.7081*n_m3(isurf)*(v_t1**4)/(nuhat*(pi**(1.5))) &
				*(dlnndpsi(isurf) + 3.3668*dlnTdpsi(isurf))
		end if

	end subroutine compute_energy_integral

! ===================================================
! Subroutine compute_P_tensor
!
! This subroutine compute the P tensor needed for
! the adjoint equation. The energy integral factor
! is not included. 
!
! ===================================================
	subroutine compute_P_tensor(isurf)

		use input_mod, only: nsurf, nlambda, nalpha, nwell
		use grids_mod, only: lambdas, alphas, dlambda
		use extrema_mod, only: min_B, max_B
		use splines_mod, only: zetas_spline, nzeta_spline
		use geometry_mod, only: compute_B, iota

		implicit none

		integer, intent(in) :: isurf
		integer :: iclass, ialpha, ilambda, izeta
		real(dp) :: lambda, lambda_scaled, theta, BB, radicand, inv_radicand

		do ilambda=1,nlambda
			lambda_scaled = lambdas(ilambda)
			lambda = 1.0/max_B(isurf) + &
				(1.0/min_B(isurf) - 1.0/max_B(isurf)) * lambda_scaled
			do ialpha=1,nalpha
				do izeta = 1,nzeta_spline
					theta = alphas(ialpha) + iota(isurf)*zetas_spline(izeta)
					BB = compute_B(isurf, theta, zetas_spline(izeta))
					radicand = max(1-lambda*BB,0.0)
					inv_radicand = max(1/(1-lambda*BB),0.0)
					do iclass=1,nclass(isurf,ilambda,ialpha)
								P_tensor_bb_before_integral(isurf,ialpha,izeta,ilambda,iclass) = &
									-2*((I_bounce_integral(isurf,ilambda,ialpha,iclass)**(-2)) &
									*(dKdalpha(isurf,ilambda,ialpha,iclass)**2)*BB*sqrt(inv_radicand)/2.0 &
									+ (-2*dIdalpha(isurf,ilambda,ialpha,iclass) &
									*dKdalpha(isurf,ilambda,ialpha,iclass) &
									*(I_bounce_integral(isurf,ilambda,ialpha,iclass)**(-2)) &
									+ 2*d2Kdalpha2(isurf,ilambda,ialpha,iclass) &
										*(I_bounce_integral(isurf,ilambda,ialpha,iclass)**(-1))) &
									*(3*BB*sqrt(radicand)/2))
								P_tensor_I_before_integral(isurf,ialpha,izeta,ilambda,iclass) = &
									2*((I_bounce_integral(isurf,ilambda,ialpha,iclass)**(-2)) &
									*(dKdalpha(isurf,ilambda,ialpha,iclass)**2)*(BB*sqrt(inv_radicand)/2.0 &
									+ sqrt(radicand)/lambda) + (-2*dIdalpha(isurf,ilambda,ialpha,iclass) &
									*dKdalpha(isurf,ilambda,ialpha,iclass) &
									*(I_bounce_integral(isurf,ilambda,ialpha,iclass)**(-2)) &
									+ 2*d2Kdalpha2(isurf,ilambda,ialpha,iclass) &
									*(I_bounce_integral(isurf,ilambda,ialpha,iclass)**(-1))) &
									*(3*BB*sqrt(radicand)/2.0 + sqrt(radicand)**3/lambda))
					end do ! iclass
				end do ! izeta
			end do ! ialpha
		end do ! ilambda

		do ialpha=1,nalpha
			do izeta=1,nzeta_spline
				! Trapezoid rule for lambda integration
				P_tensor_bb_before_integral(isurf,ialpha,izeta,nlambda,:) = &
					0.5*P_tensor_bb_before_integral(isurf,ialpha,izeta,nlambda,:)
				P_tensor_bb_before_integral(isurf,ialpha,izeta,1,:) = &
					0.5*P_tensor_bb_before_integral(isurf,ialpha,izeta,1,:)
				P_tensor_I_before_integral(isurf,ialpha,izeta,nlambda,:) = &
					0.5*P_tensor_I_before_integral(isurf,ialpha,izeta,nlambda,:)
				P_tensor_I_before_integral(isurf,ialpha,izeta,1,:) = &
					0.5*P_tensor_I_before_integral(isurf,ialpha,izeta,1,:)
				! Perform integration over lambda and sum over particle class
				P_tensor_bb(isurf,ialpha,izeta) = sum(P_tensor_bb_before_integral(isurf,ialpha,izeta,:,:))
				P_tensor_I(isurf,ialpha,izeta) = sum(P_tensor_I_before_integral(isurf,ialpha,izeta,:,:))
				P_tensor_bb(isurf,ialpha,izeta) = P_tensor_bb(isurf,ialpha,izeta)*dlambda(isurf)
				P_tensor_I(isurf,ialpha,izeta) = P_tensor_I(isurf,ialpha,izeta)*dlambda(isurf)
			end do
		end do

	end subroutine compute_P_tensor

end module diagnostics_mod
