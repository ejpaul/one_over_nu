! ===================================================
! subroutine write_output
!
! This subroutine writes to a netcdf output file
! using the ezcdf interface.
!
! ===================================================
subroutine write_output(total_time)

	use stel_kinds
	use ezcdf
	use grids_mod, only: alphas, lambdas
	use extrema_mod, only: min_B, max_B, B, thetas, zetas
	use geometry_mod, only: xm, xn, s_surf, ns, nmodes, nfp, iota, Boozer_I, Boozer_G, bmnc, &
			B_pest_mnc, mnmax_transform, xm_transform, xn_transform, B_pest_mnc, B_dot_grad_zeta_pest_mnc
	use input_mod
	use diagnostics_mod

	implicit none

	real(dp), intent(in) :: total_time
	integer :: ierr, ncid
	character(len=200) :: output_filename

	! Scalars:
	character(len=*), parameter :: &
		vn_nfp = "nfp", &
		vn_nlambda = "nlambda", &
		vn_nalpha = "nalpha", &
		vn_ntheta = "ntheta", &
		vn_nzeta = "nzeta", &
		vn_nintegral = "nintegral", &
		vn_Delta_zeta = "Delta_zeta", &
		vn_max_search_in_zeta = "max_search_in_zeta", &
		vn_ns = "ns", &
		vn_nmodes = "nmodes", &
		vn_nwell = "nwell", &
		vn_niter_newton = "niter_newton", &
		vn_tol_newton = "tol_newton", &
		vn_root_search_tolerance = "root_search_tolerance", &
		vn_Niter_root = "Niter_root", &
		vn_total_time = "total_time", &
		vn_output_particle_flux  = "output_particle_flux", &
		vn_collision_species_option = "collision_species_option"

	! Arrays with dimension 1
	character(len=*), parameter :: &
		vn_thetas = "thetas", &
		vn_zetas = "zetas", &
		vn_xm = "xm", &
		vn_xn = "xn", &
		vn_lambdas = "lambdas", &
		vn_alphas = "alphas", &
		vn_boozmn_filename = "boozmn_filename", &
		vn_one_over_nu_metric = "one_over_nu_metric", &
		vn_min_B = "min_B", &
		vn_max_B = "max_B", &
		vn_iota = "iota", &
		vn_Boozer_I = "Boozer_I", &
		vn_Boozer_G = "Boozer_G", &
		vn_s_wish = "s_wish", &
		vn_s_surf = "s_surf", &
		vn_particleFlux = "particleFlux", &
		vn_energy_integral = "energy_integral", &
		vn_wout_filename = "wout_filename", &
		vn_xm_transform = "xm_transform", &
		vn_xn_transform = "xn_transform"

	! Arrays with dimension 2
	character(len=*), parameter :: &
		vn_bmnc = "bmnc", &
		vn_b_pest_mnc = "B_pest_mnc", &
		vn_B_dot_grad_zeta_pest_mnc = "B_dot_grad_zeta_pest_mnc"

	! Arrays with dimension 3
	character(len=*), parameter :: &
		vn_B = "B", &
		vn_nclass = "nclass", &
		vn_P_tensor_bb = "P_tensor_bb", &
		vn_P_tensor_I = "P_tensor_I"

	! Arrays with dimension 4 (summed over nwell dimension)
	character(len=*), parameter :: &
		vn_J_invariant = "J_invariant", &
		vn_dKdalpha = "dKdalpha", &
		vn_I_bounce_integral = "I_bounce_integral", &
		vn_one_over_nu_metric_before_integral = "one_over_nu_metric_before_integral", &
		vn_K_bounce_integral = "K_bounce_integral"

	! Arrays with dimension 1
	character(len=*), parameter, dimension(1) :: &
		ntheta_dim = (/'ntheta'/), &
		nzeta_dim = (/'nzeta'/), &
		nmodes_dim = (/'nmodes'/), &
		nlambda_dim = (/'nlambda'/), &
		nalpha_dim = (/'nalpha'/), &
		string_length_dim = (/'string_length'/), &
		nsurf_dim = (/'nsurf'/), &
		mnmax_transform_dim = (/'mnmax_transform'/)

	! Arrays with dimension 2
	character(len=*), parameter, dimension(2) :: &
		nsurf_nmodes_dim = (/character(len=50) :: 'nsurf','nmodes' /), &
		nsurf_mnmax_transform_dim = (/ character(len=50) :: &
			'nsurf','mnmax_transform' /)

	! Arrays with dimension 3
	character(len=*), parameter, dimension(3) :: &
		nsurf_ntheta_nzeta_dim = (/ character(len=50) :: 'nsurf','ntheta','nzeta'/), &
		nsurf_nlambda_nalpha_dim = (/ character(len=50) :: 'nsurf','nlambda', 'nalpha'/), &
		nsurf_nalpha_nzeta_spline_dim = (/ character(len=50) :: &
			'nsurf','nalpha','nzeta_spline'/)

	output_filename = "one_over_nu_out" // trim(inputFilename(15:)) // ".nc"

	call cdf_open(ncid,output_filename,'w',ierr)
	IF (ierr .ne. 0) then
		print *,"Error opening output file ",output_filename
		stop
	end IF

	! Scalars
	call cdf_define(ncid, vn_nfp, nfp)
	call cdf_define(ncid, vn_nlambda, nlambda)
	call cdf_define(ncid, vn_nalpha, nalpha)
	call cdf_define(ncid, vn_ntheta, ntheta)
	call cdf_define(ncid, vn_nzeta, nzeta)
	call cdf_define(ncid, vn_nintegral, nintegral)
	call cdf_define(ncid, vn_Delta_zeta, Delta_zeta)
	call cdf_define(ncid, vn_max_search_in_zeta, max_search_in_zeta)
	call cdf_define(ncid, vn_ns, ns)
	call cdf_define(ncid, vn_nmodes, nmodes)
	call cdf_define(ncid, vn_nwell, nwell)
	call cdf_define(ncid, vn_niter_newton, niter_newton)
	call cdf_define(ncid, vn_tol_newton, tol_newton)
	call cdf_define(ncid, vn_root_search_tolerance, root_search_tolerance)
	call cdf_define(ncid, vn_niter_root, niter_root)
	call cdf_define(ncid, vn_total_time, total_time)
	call cdf_define(ncid, vn_output_particle_flux, output_particle_flux)
	call cdf_define(ncid, vn_collision_species_option, collision_species_option)

	! 1 dimension
	call cdf_define(ncid,vn_thetas,thetas,dimname=ntheta_dim)
	call cdf_define(ncid,vn_zetas,zetas,dimname=nzeta_dim)
	call cdf_define(ncid,vn_s_wish,s_wish(1:nsurf),dimname=nsurf_dim)
	if (geometry_option == 1) then
		call cdf_define(ncid,vn_xm,xm,dimname=nmodes_dim)
		call cdf_define(ncid,vn_xn,xn,dimname=nmodes_dim)
	else
		call cdf_define(ncid,vn_xm_transform,xm_transform,dimname=mnmax_transform_dim)
		call cdf_define(ncid,vn_xn_transform,xn_transform,dimname=mnmax_transform_dim)
	end if
	call cdf_define(ncid,vn_lambdas,lambdas,dimname=nlambda_dim)
	call cdf_define(ncid,vn_alphas,alphas,dimname=nalpha_dim)
	if (geometry_option==1) then
		call cdf_define(ncid,vn_boozmn_filename,boozmn_filename,dimname=string_length_dim)
	else
		call cdf_define(ncid,vn_wout_filename,wout_filename,dimname=string_length_dim)
	end if
	call cdf_define(ncid, vn_s_surf, s_surf, dimname=nsurf_dim)
	call cdf_define(ncid, vn_min_B, min_B, dimname=nsurf_dim)
	call cdf_define(ncid, vn_max_B, max_B, dimname=nsurf_dim)
	call cdf_define(ncid, vn_one_over_nu_metric, one_over_nu_metric, dimname=nsurf_dim)
	call cdf_define(ncid, vn_iota, iota, dimname=nsurf_dim)
	if (geometry_option == 1) then
		call cdf_define(ncid, vn_Boozer_G, Boozer_G, dimname=nsurf_dim)
		call cdf_define(ncid, vn_Boozer_I, Boozer_I, dimname=nsurf_dim)
	end if
	if (output_particle_flux) then
		call cdf_define(ncid, vn_particleFlux, particleFlux, dimname=nsurf_dim)
	end if
	if (output_p_tensor .or. output_particle_flux) then
		call cdf_define(ncid, vn_energy_integral,energy_integral, dimname=nsurf_dim)
	end if

	! 2 dimension
	if (geometry_option==1) then
		call cdf_define(ncid,vn_bmnc,bmnc,dimname=nsurf_nmodes_dim)
	else
		call cdf_define(ncid,vn_B_pest_mnc,B_pest_mnc,dimname=nsurf_mnmax_transform_dim)
		call cdf_define(ncid,vn_B_dot_grad_zeta_pest_mnc,B_dot_grad_zeta_pest_mnc, &
			dimname=nsurf_mnmax_transform_dim)
	end if

	! 3 dimension
	call cdf_define(ncid,vn_B,B,dimname=nsurf_ntheta_nzeta_dim)
	call cdf_define(ncid,vn_nclass,nclass,dimname=nsurf_nlambda_nalpha_dim)
	if (output_P_tensor) then
		call cdf_define(ncid,vn_P_tensor_bb,P_tensor_bb,dimname=nsurf_nalpha_nzeta_spline_dim)
		call cdf_define(ncid,vn_P_tensor_I,P_tensor_I,dimname=nsurf_nalpha_nzeta_spline_dim)
	end if

	! 4 dimension (summed over nwell dimension)
	if (output_J) then
		call cdf_define(ncid,vn_J_invariant,sum(J_invariant,4),dimname=nsurf_nlambda_nalpha_dim)
	end if
	call cdf_define(ncid,vn_dKdalpha,sum(dKdalpha,4),dimname=nsurf_nlambda_nalpha_dim)
	call cdf_define(ncid,vn_I_bounce_integral,sum(I_bounce_integral,4),dimname=nsurf_nlambda_nalpha_dim)
	call cdf_define(ncid,vn_one_over_nu_metric_before_integral,sum(one_over_nu_metric_before_integral,4),dimname=nsurf_nlambda_nalpha_dim)

	! scalars
	call cdf_write(ncid, vn_nfp, nfp)
	call cdf_write(ncid, vn_nlambda, nlambda)
	call cdf_write(ncid, vn_nalpha, nalpha)
	call cdf_write(ncid, vn_ntheta, ntheta)
	call cdf_write(ncid, vn_nzeta, nzeta)
	call cdf_write(ncid, vn_nintegral, nintegral)
	call cdf_write(ncid, vn_Delta_zeta, Delta_zeta)
	call cdf_write(ncid, vn_max_search_in_zeta, max_search_in_zeta)
	call cdf_write(ncid, vn_ns, ns)
	if (geometry_option==1) then
		call cdf_write(ncid, vn_nmodes, nmodes)
	end if
	call cdf_write(ncid, vn_nwell, nwell)
	call cdf_write(ncid, vn_niter_newton, niter_newton)
	call cdf_write(ncid, vn_tol_newton, tol_newton)
	call cdf_write(ncid, vn_root_search_tolerance, root_search_tolerance)
	call cdf_write(ncid, vn_niter_root, niter_root)
	call cdf_write(ncid, vn_total_time, total_time)
	call cdf_write(ncid, vn_output_particle_flux, output_particle_flux)
	call cdf_write(ncid, vn_collision_species_option, collision_species_option)

	! 1 dimension
	call cdf_write(ncid, vn_s_wish, s_wish(1:nsurf))
	call cdf_write(ncid, vn_one_over_nu_metric, one_over_nu_metric)
	call cdf_write(ncid, vn_max_B, max_B)
	call cdf_write(ncid, vn_min_B, min_B)
	call cdf_write(ncid, vn_thetas,thetas)
	call cdf_write(ncid, vn_zetas,zetas)
	call cdf_write(ncid, vn_iota, iota)
	if (geometry_option==1) then
		call cdf_write(ncid, vn_Boozer_G, Boozer_G)
		call cdf_write(ncid, vn_Boozer_I, Boozer_I)
		call cdf_write(ncid,vn_xm,xm)
		call cdf_write(ncid,vn_xn,xn)
	else
		call cdf_write(ncid, vn_B_pest_mnc, B_pest_mnc)
		call cdf_write(ncid, vn_B_dot_grad_zeta_pest_mnc, B_dot_grad_zeta_pest_mnc)
		call cdf_write(ncid, vn_xm_transform, xm_transform)
		call cdf_write(ncid, vn_xn_transform, xn_transform)
	end if
	call cdf_write(ncid,vn_thetas,thetas)
	call cdf_write(ncid,vn_zetas,zetas)
	call cdf_write(ncid,vn_lambdas,lambdas)
	call cdf_write(ncid,vn_alphas,alphas)
	if (geometry_option == 1) then
		call cdf_write(ncid,vn_boozmn_filename,trim(boozmn_filename))
	else
		call cdf_write(ncid,vn_wout_filename,trim(wout_filename))
	end if
	call cdf_write(ncid,vn_s_surf,s_surf)
	if (output_particle_flux) then
		call cdf_write(ncid,vn_particleFlux,particleFlux)
	end if
	if (output_particle_flux .or. output_p_tensor) then
		call cdf_write(ncid,vn_energy_integral,energy_integral)
	end if

	! 2 dimension
	if (geometry_option == 1) then
		call cdf_write(ncid,vn_bmnc,bmnc)
	end if

	! 3 dimension
	call cdf_write(ncid,vn_B,B)
	call cdf_write(ncid,vn_nclass,nclass)
	if (output_P_tensor) then
		call cdf_write(ncid,vn_P_tensor_bb,P_tensor_bb)
		call cdf_write(ncid,vn_P_tensor_I,P_tensor_I)
	end if

	! 4 dimension (summed over nwell)
	if (output_J) then
		call cdf_write(ncid,vn_J_invariant,sum(J_invariant,4))
	end if 
	call cdf_write(ncid,vn_dKdalpha,sum(dKdalpha,4))
	call cdf_write(ncid,vn_I_bounce_integral,sum(I_bounce_integral,4))
	call cdf_write(ncid,vn_one_over_nu_metric_before_integral,sum(one_over_nu_metric_before_integral,4))

	call cdf_close(ncid)


end subroutine write_output
