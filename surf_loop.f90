 ! ===================================================
 ! Subroutine surf_loop
 !
 ! This subroutine loops over the desired magnetic
 ! surfaces and lambda grid to compute the bounce
 ! bounce integrals. The 1/nu metric is computed
 ! by summing over particle classes and integrating
 ! over alpha and lambda.
 !
 ! ===================================================
subroutine surf_loop()

	use stel_constants
	use flintegrate_mod, only: flintegrate
	use omp_lib
	use diagnostics_mod, only: compute_diagnostics, init_diagnostics
	use input_mod, only: verbose, nlambda, nsurf

	implicit none

	integer :: isurf, ilambda

	call init_diagnostics

	do isurf = 1,nsurf
		do ilambda = 1,nlambda
			if (verbose) then
				print *,"ilambda: ", ilambda
			end if
			call flintegrate(isurf,ilambda)
		end do
		call compute_diagnostics(isurf)
	end do ! isurf

end subroutine surf_loop



