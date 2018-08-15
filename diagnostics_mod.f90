module diagnostics_mod

	use stel_kinds

	! Bounce integrals
	real(dp), dimension(:,:,:,:), allocatable :: J_invariant, dKdalpha, I_bounce_integral, one_over_nu_metric_before_integral, K_bounce_integral, nemov_metric_before_integral, H_bounce_integral
	! Number of classes for each grid point
	integer, dimension(:,:,:), allocatable :: nclass
	! Output metrics
	real(dp), dimension(:), allocatable :: one_over_nu_metric, nemov_metric

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
		use input_mod, only: nlambda, nalpha

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

	end subroutine compute_diagnostics

end module diagnostics_mod
