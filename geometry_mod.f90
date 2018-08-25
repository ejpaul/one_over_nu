module geometry_mod

	use constants_mod
	use stel_kinds
	use input_mod, only: boozmn_filename, nsurf, wout_filename, geometry_option
	use read_wout_mod, only: xm_vmec => xm, xn_vmec => xn

	! Boozmn geometry items
	integer :: ns, nfp, nmodes
	real(dp), dimension(:), allocatable :: iota, Boozer_I, Boozer_G, xm, xn
	real(dp), dimension(:,:), allocatable :: bmnc
	! Value of s for nearest surface in boozmn_filename (used for calculation)
	real(dp), dimension(:), allocatable :: s_surf
	! PEST coordinate items
	real(dp), dimension(:,:), allocatable :: B_pest_mnc, B_dot_grad_zeta_pest_mnc, &
		Bsubu_pest_mnc, Bsubv_pest_mnc
	integer, dimension(:), allocatable :: xm_transform, xn_transform
	integer :: mnmax_transform

	! Needed for call to FZERO
	real(dp) :: theta_pest_target, this_zeta
	real(dp), dimension(:), allocatable :: lmns_slice
	private :: theta_pest_target, this_zeta, lmns_slice

	contains

! ===================================================
! function compute_geometry
!
! The function computes B, dBdtheta, dBdzeta for given
! values of thet and zeta using the bmnc's.
!
!	Inputs:
! isurf			index in ssurf corresponding to calculation
! theta			value of theta for calculation
! zeta			value of zeta for calculation
!
!	Outputs:
!	compute_geometry	array of length geometry_length
!										This array can be indexed into
!										using parameters defined
!										in constants_mod.
!
! ===================================================
	function compute_geometry(isurf,theta,zeta)

		use stel_constants

		real(dp) :: theta,zeta
		real(dp), dimension(geometry_length) :: compute_geometry
		integer :: imn
		real(dp), dimension(:), allocatable :: cosangle, sinangle
		integer, intent(in) :: isurf
		real(dp) :: BB, dBBdtheta, sign_BBdotgradzeta, BBdotgradzeta

		if (geometry_option == 1) then
			allocate(cosangle(nmodes))
			allocate(sinangle(nmodes))
			cosangle = cos(xm*theta-xn*zeta)
			sinangle = sin(xm*theta-xn*zeta)
			BB = dot_product(bmnc(isurf,:),cosangle)
			compute_geometry(B_index) = BB
			dBBdtheta = -dot_product(bmnc(isurf,:),xm*sinangle)
			compute_geometry(dBdtheta_index) = dBBdtheta
			compute_geometry(Bdotgradzeta_index) = abs((BB**2) &
				/(Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf)))
			sign_BBdotgradzeta = sign(one,Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf))
			compute_geometry(dBdotgradzetadtheta_index) = (2*BB*dBBdtheta &
				/(Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf)))*sign_BBdotgradzeta
		else
			allocate(cosangle(mnmax_transform))
			allocate(sinangle(mnmax_transform))
			cosangle = cos(xm_transform*theta-xn_transform*zeta)
			sinangle = sin(xm_transform*theta-xn_transform*zeta)
			compute_geometry(B_index) = dot_product(B_pest_mnc(isurf,:),cosangle)
			compute_geometry(dBdtheta_index) = -dot_product(B_pest_mnc(isurf,:),xm_transform*sinangle)
			BBdotgradzeta = dot_product(B_dot_grad_zeta_pest_mnc(isurf,:), &
				cosangle)
			compute_geometry(Bdotgradzeta_index) = abs(BBdotgradzeta)
			sign_BBdotgradzeta = sign(one,BBdotgradzeta)
			compute_geometry(dBdotgradzetadtheta_index) = -dot_product(B_dot_grad_zeta_pest_mnc(isurf,:) &
				*xm_transform, sinangle)*sign_BBdotgradzeta
		end if

	end function

	real(dp) function compute_B(isurf,theta,zeta)

		real(dp) :: theta,zeta
		integer :: isurf

		if (geometry_option == 1) then
			compute_B = dot_product(bmnc(isurf,:),cos(xm*theta-xn*zeta))
		else
			compute_B = dot_product(B_pest_mnc(isurf,:),cos(xm_transform*theta-xn_transform*zeta))
		end if

	end function compute_B

	real(dp) function compute_d2Bdtheta2(isurf,theta,zeta)

		real(dp) :: theta, zeta
		integer :: isurf

		if (geometry_option==1) then
			compute_d2Bdtheta2 = -dot_product(bmnc(isurf,:)*xm*xm,cos(xm*theta-xn*zeta))
		else
			compute_d2Bdtheta2 = -dot_product(B_pest_mnc(isurf,:)*xm_transform*xm_transform,&
				cos(xm_transform*theta-xn_transform*zeta))
		end if

	end function compute_d2Bdtheta2

	function compute_dBdzeta(isurf,theta,zeta)

		real(dp), intent(in) :: theta, zeta
		integer, intent(in) :: isurf
		real(dp) :: compute_dBdzeta

		if (geometry_option==1) then
			compute_dBdzeta = dot_product(bmnc(isurf,:)*xn,sin(xm*theta-xn*zeta))
		else
			compute_dBdzeta = dot_product(B_pest_mnc(isurf,:)*xn_transform,&
				sin(xm_transform*theta-xn_transform*zeta))
		end if

	end function compute_dBdzeta

	real(dp) function compute_d2Bdzeta2(isurf,theta,zeta)

		real(dp) :: theta, zeta
		integer :: isurf

		if (geometry_option==1) then
			compute_d2Bdzeta2 = -dot_product(bmnc(isurf,:)*xn*xn,cos(xm*theta-xn*zeta))
		else
			compute_d2Bdzeta2 = -dot_product(B_pest_mnc(isurf,:)*xn_transform*xn_transform,&
				cos(xm_transform*theta-xn_transform*zeta))
		end if

	end function compute_d2Bdzeta2

	real(dp) function compute_d2Bdthetadzeta(isurf,theta,zeta)

		real(dp) :: theta, zeta
		integer :: isurf

		if (geometry_option==1) then
			compute_d2Bdthetadzeta = dot_product(bmnc(isurf,:)*xm*xn,cos(xm*theta-xn*zeta))
		else
			compute_d2Bdthetadzeta = dot_product(B_pest_mnc(isurf,:)*xm_transform*xn_transform,&
				cos(xm_transform*theta-xn_transform*zeta))
		end if

	end function compute_d2Bdthetadzeta

	function compute_d2Bdotgradzetadtheta2(isurf,theta,zeta)

		use stel_constants

		real(dp), intent(in) :: theta, zeta
		integer, intent(in) :: isurf
		real(dp) :: compute_d2Bdotgradzetadtheta2
		real(dp) :: BB, dBBdtheta, d2BBdtheta2, sign_Bdotgradzeta, BBdotgradzeta

		if (geometry_option==1) then
			BB = dot_product(bmnc(isurf,:),cos(xm*theta-xn*zeta))
			dBBdtheta = -dot_product(bmnc(isurf,:)*xm,sin(xm*theta-xn*zeta))
			d2BBdtheta2 = -dot_product(bmnc(isurf,:)*xm*xm,cos(xm*theta-xn*zeta))
			sign_Bdotgradzeta = sign(one,Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf))
			compute_d2Bdotgradzetadtheta2 = 2*(BB*d2BBdtheta2 + dBBdtheta**2) &
				/(Boozer_G(isurf) + iota(isurf)*Boozer_I(isurf))*sign_Bdotgradzeta
		else
			BBdotgradzeta = dot_product(B_dot_grad_zeta_pest_mnc(isurf,:),&
				cos(xm_transform*theta-xn_transform*zeta))
			sign_Bdotgradzeta = sign(one,BBdotgradzeta)
			compute_d2Bdotgradzetadtheta2 = dot_product(B_dot_grad_zeta_pest_mnc(isurf,:)*xm_transform*xm_transform,&
				cos(xm_transform*theta-xn_transform*zeta))*sign_Bdotgradzeta
		end if

	end function compute_d2Bdotgradzetadtheta2

! ===================================================
! subroutine init_geometry_boozer
!
!	This subroutine reads the boozmn input file
! to get the required geometric quantities on surfaces
! closest to the on requested surfaces in s_wish. This
! uses the read_boozer_mod from LIBSTELL.
!
! ===================================================
	subroutine init_geometry_boozer

		use read_boozer_mod
		use stel_constants
		use input_mod, only: verbose, s_wish

		implicit none

		integer :: ierr, iopen, is, ns_avail, s_index, s_avail, jsize, ilist, i
		integer :: itheta, izeta, imn, ialpha, ilambda
		real(dp), dimension(:), allocatable :: s_half, s_half_available
		integer, dimension(:), allocatable :: jlist
		real(dp) :: ds, ds_curr, s, angle
		integer :: isurf

		call read_boozer_file(boozmn_filename, ierr, iopen)
		if (iopen .ne. 0) stop 'Error opening boozmn file'
		if (ierr .ne. 0) stop 'Error reading boozmn file'
		if (verbose) print *,"  Successfully read boozmn data from ",trim(boozmn_filename)

		allocate(bmnc(nsurf,mnboz_b))
		allocate(Boozer_I(nsurf))
		allocate(Boozer_G(nsurf))
		allocate(iota(nsurf))
		allocate(s_surf(nsurf))

		nfp = nfp_b
		ns = ns_b

		! Find available s in boozmn file
		ds = 1.0/(ns-1)

		allocate(s_half(ns-1))
		s_half = 0
		do is=2,ns
			s_half(is-1) = (is-1)*ds - 0.5*ds
		end do
		jsize = count(idx_b ==1)
		i = 1
		allocate(jlist(jsize))
		allocate(s_half_available(jsize))
		do ilist = 1,ns
			if (idx_b(ilist) .le. 0) cycle
			jlist(i) = ilist
			s_half_available(i) = s_half(ilist-1)
			i = i + 1
		end do

		allocate(xm(mnboz_b))
		allocate(xn(mnboz_b))
		xm = ixm_b
		xn = ixn_b
		nmodes = mnboz_b

		do isurf=1,nsurf
			! Find closest available s in boozmn file
			ds_curr = 1.0
			s_index = 1
			do is=1,jsize
				if (abs(s_half_available(is)-s_wish(isurf))<ds_curr) then
					s = s_half_available(is)
					ds_curr = abs(s_half_available(is)-s_wish(isurf))
					s_index = is
				end if
			end do
			s_surf(isurf) = s_half_available(s_index)
			Boozer_I(isurf) = buco_b(jlist(s_index))
			Boozer_G(isurf) = bvco_b(jlist(s_index))
			iota(isurf) = iota_b(jlist(s_index))
			bmnc(isurf,:) = bmnc_b(:,jlist(s_index))

		end do ! isurf

	end subroutine init_geometry_boozer

	subroutine init_geometry_vmec

		use read_wout_mod, only: nfp_vmec => nfp,  &
			ns_vmec => ns, read_wout_file, iota_vmec => iotas, bmnc_vmec => bmnc, &
			mnmax_nyq, gmnc_vmec => gmnc, lmns_vmec => lmns, bsubumnc, bsubvmnc, &
			xm_nyq_vmec => xm_nyq, xn_nyq_vmec => xn_nyq, mpol_vmec => mpol, &
			ntor_vmec => ntor, mnmax, phi_vmec => phi, xm_vmec => xm, xn_vmec => xn
		use input_mod, only: verbose, s_wish, mpol_transform_refinement, &
			ntor_transform_refinement, transform_relerr, transform_abserr
		use stel_constants

		implicit none

		integer :: ierr, iopen, isurf, is, s_index, itheta, izeta, mpol_transform, &
			ntor_transform, ntheta_transform, nzeta_transform
		real(dp), dimension(:), allocatable :: s_full, s_half, bmnc_slice, gmnc_slice, &
			 bsubumnc_slice, bsubvmnc_slice, thetas_transform, zetas_transform
		real(dp) :: theta_min, theta_max, ds
		real(dp), dimension(:,:), allocatable :: theta_vmec_on_theta_pest_grid, &
			B_on_theta_pest_grid, Bsubu_on_theta_pest_grid, Bsubv_on_theta_pest_grid, &
			Bsubs_on_theta_pest_grid, sqrt_g_vmec_on_theta_pest_grid, &
			d_lambda_d_theta_vmec_on_theta_pest_grid, lambda_on_theta_vmec_grid, &
			B_dot_grad_zeta_on_theta_pest_grid, B_on_theta_pest_grid_reconstructed, &
			Bsubu_on_theta_pest_grid_reconstructed, Bsubv_on_theta_pest_grid_reconstructed, &
			B_dot_grad_zeta_on_theta_pest_grid_reconstructed
		integer :: fzero_flag, im, jn, index, imn
		real(dp) :: dpsids, factor, factor2
		real(dp), dimension(:,:), allocatable :: cosangle, sinangle, thetas_transform_2d, &
			zetas_transform_2d, angle

		allocate(iota(nsurf))
		allocate(s_surf(nsurf))

		call read_wout_file(wout_filename, ierr, iopen)
		if (iopen .ne. 0) stop 'Error opening wout file'
		if (ierr .ne. 0) stop 'Error reading wout file'
		if (verbose) print *,"  Successfully read wout data from ",trim(wout_filename)

		nfp = nfp_vmec
		ns = ns_vmec

		allocate(s_full(ns))
		allocate(s_half(ns))
		do is =1,ns
			s_full(is) = 1.0*(is-1.0)/(ns-1)
		end do
		ds = s_full(2) - s_full(1)
		s_half = s_full(2:ns) - 0.5*ds
		dpsids = phi_vmec(ns)/(twopi)
		allocate(bmnc_slice(mnmax_nyq))
		allocate(gmnc_slice(mnmax_nyq))
		allocate(lmns_slice(mnmax_nyq))
		allocate(bsubumnc_slice(mnmax_nyq))
		allocate(bsubvmnc_slice(mnmax_nyq))

		! Construct theta and zeta grids for coordinate transformation
		mpol_transform = mpol_vmec*mpol_transform_refinement
		ntor_transform = ntor_vmec*ntor_transform_refinement
		ntheta_transform = mpol_transform*2
		nzeta_transform = ntor_transform*2
		mnmax_transform = (2*ntor_transform+1)*mpol_transform + ntor_transform+1
		allocate(xm_transform(mnmax_transform))
		allocate(xn_transform(mnmax_transform))
		xm_transform = zero
		xn_transform = zero
		! Initialize m = 0 modes
		do jn = 1,ntor_transform+1
			xn_transform(jn) = (jn-1)*nfp
		end do
		! Initialize m > 0 modes
		index = ntor_transform+2
		do im = 1,mpol_transform
			do jn = -ntor_transform*nfp,ntor_transform*nfp,nfp
				xm_transform(index) = im
				xn_transform(index) = jn
				index = index + 1
			end do
		end do

		allocate(B_pest_mnc(nsurf,mnmax_transform))
		allocate(B_dot_grad_zeta_pest_mnc(nsurf,mnmax_transform))
		allocate(Bsubu_pest_mnc(nsurf,mnmax_transform))
		allocate(Bsubv_pest_mnc(nsurf,mnmax_transform))

		allocate(thetas_transform(ntheta_transform))
		allocate(zetas_transform(nzeta_transform))
		allocate(thetas_transform_2d(ntheta_transform,nzeta_transform))
		allocate(zetas_transform_2d(ntheta_transform,nzeta_transform))
		do itheta=1,ntheta_transform
			thetas_transform(itheta) = twopi*(itheta-1.0)/ntheta_transform
			thetas_transform_2d(itheta,:) = thetas_transform(itheta)
		end do
		do izeta=1,nzeta_transform
			zetas_transform(izeta) = (twopi/nfp)*(izeta-1.0)/nzeta_transform
			zetas_transform_2d(:,izeta) = zetas_transform(izeta)
		end do

		allocate(angle(ntheta_transform,nzeta_transform))
		allocate(cosangle(ntheta_transform,nzeta_transform))
		allocate(sinangle(ntheta_transform,nzeta_transform))

		allocate(theta_vmec_on_theta_pest_grid(ntheta_transform,nzeta_transform))

		do isurf=1,nsurf
			! Find closest available half-grid surface to desired value of s
			s_index = minloc(abs(s_half-s_wish(isurf)),1)
			s_surf(isurf) = s_half(s_index)
			! Add 1 to account for the 0 at the beginning of half-mesh quantities.
			s_index = s_index + 1
			! Half mesh quantities
			iota(isurf) =	iota_vmec(s_index)
			bmnc_slice = bmnc_vmec(:,s_index)
			gmnc_slice = gmnc_vmec(:,s_index)
			lmns_slice = lmns_vmec(:,s_index)
			bsubumnc_slice = bsubumnc(:,s_index)
			bsubvmnc_slice = bsubvmnc(:,s_index)

			! Perform root solve for theta_vmec_on_theta_pest_grid
			theta_vmec_on_theta_pest_grid = zero
			do itheta=1,ntheta_transform
				theta_pest_target = thetas_transform(itheta)
				do izeta = 1, nzeta_transform
					this_zeta = zetas_transform(izeta)
					! Initial point for root solve
					theta_min = theta_pest_target - 0.3
					theta_max = theta_pest_target + 0.3
					call fzero(fzero_residual, theta_min, theta_max, theta_pest_target, &
						transform_relerr, transform_abserr, fzero_flag)
					theta_vmec_on_theta_pest_grid(itheta,izeta) = theta_min
					 if (fzero_flag == 4) then
							stop "ERROR: fzero returned error 4: no sign change in residual"
					 else if (fzero_flag > 2) then
							print *,"WARNING in irp: fzero returned an error code:",fzero_flag
					 end if
				end do
			end do

			! Evaluate all quantities using uniformly spaced grid in theta_pest
			allocate(B_on_theta_pest_grid(ntheta_transform,nzeta_transform))
			allocate(Bsubu_on_theta_pest_grid(ntheta_transform,nzeta_transform))
			allocate(Bsubv_on_theta_pest_grid(ntheta_transform,nzeta_transform))
			allocate(sqrt_g_vmec_on_theta_pest_grid(ntheta_transform,nzeta_transform))
			allocate(d_lambda_d_theta_vmec_on_theta_pest_grid(ntheta_transform,nzeta_transform))
			allocate(lambda_on_theta_vmec_grid(ntheta_transform,nzeta_transform))
			allocate(B_dot_grad_zeta_on_theta_pest_grid(ntheta_transform,nzeta_transform))
			B_on_theta_pest_grid = zero
			Bsubu_on_theta_pest_grid = zero
			Bsubv_on_theta_pest_grid = zero
			sqrt_g_vmec_on_theta_pest_grid = zero
			d_lambda_d_theta_vmec_on_theta_pest_grid = zero
			lambda_on_theta_vmec_grid = zero
			B_dot_grad_zeta_on_theta_pest_grid = zero

			! Reconstructed quantities
			allocate(B_on_theta_pest_grid_reconstructed(ntheta_transform,nzeta_transform))
			allocate(Bsubu_on_theta_pest_grid_reconstructed(ntheta_transform,nzeta_transform))
			allocate(Bsubv_on_theta_pest_grid_reconstructed(ntheta_transform,nzeta_transform))
			allocate(B_dot_grad_zeta_on_theta_pest_grid_reconstructed(ntheta_transform,nzeta_transform))
			B_on_theta_pest_grid_reconstructed = zero
			Bsubu_on_theta_pest_grid_reconstructed = zero
			Bsubv_on_theta_pest_grid_reconstructed = zero
			B_dot_grad_zeta_on_theta_pest_grid_reconstructed = zero

			! Compute quantities on theta_pest grid
			do itheta=1,ntheta_transform
				do izeta=1,nzeta_transform
				! Quantities that use non-nyquist xm and xn modes
					lambda_on_theta_vmec_grid(itheta,izeta) = lambda_on_theta_vmec_grid(itheta,izeta) &
						+ dot_product(lmns_slice,sin(xm_vmec*thetas_transform(itheta)&
						-xn_vmec*zetas_transform(izeta)))
					d_lambda_d_theta_vmec_on_theta_pest_grid(itheta,izeta) = &
						d_lambda_d_theta_vmec_on_theta_pest_grid(itheta,izeta) &
						+ dot_product(lmns_slice*xm_vmec,cos(xm_vmec*theta_vmec_on_theta_pest_grid(itheta,izeta) &
						- xn_vmec*zetas_transform(izeta)))
				! Quantities that use nyquist xm and xn modes
					B_on_theta_pest_grid(itheta,izeta) = B_on_theta_pest_grid(itheta,izeta) &
						+ dot_product(bmnc_slice,cos(xm_nyq_vmec*theta_vmec_on_theta_pest_grid(itheta,izeta) &
						- xn_nyq_vmec*zetas_transform(izeta)))
					Bsubu_on_theta_pest_grid(itheta,izeta) = Bsubu_on_theta_pest_grid(itheta,izeta) &
						+ dot_product(bsubumnc_slice,cos(xm_nyq_vmec*theta_vmec_on_theta_pest_grid(itheta,izeta) &
						- xn_nyq_vmec*zetas_transform(izeta)))
					Bsubv_on_theta_pest_grid(itheta,izeta) = Bsubv_on_theta_pest_grid(itheta,izeta) &
						+ dot_product(bsubvmnc_slice,cos(xm_nyq_vmec*theta_vmec_on_theta_pest_grid(itheta,izeta) &
						- xn_nyq_vmec*zetas_transform(izeta)))
					sqrt_g_vmec_on_theta_pest_grid(itheta,izeta) = sqrt_g_vmec_on_theta_pest_grid(itheta,izeta) &
						+ dot_product(gmnc_slice,cos(xm_nyq_vmec*theta_vmec_on_theta_pest_grid(itheta,izeta) &
						- xn_nyq_vmec*zetas_transform(izeta)))
				end do
			end do
			B_dot_grad_zeta_on_theta_pest_grid = dpsids*(1.0 + d_lambda_d_theta_vmec_on_theta_pest_grid) &
				/ sqrt_g_vmec_on_theta_pest_grid
			! Fourier transform
			B_pest_mnc(isurf,1) = sum(B_on_theta_pest_grid)/(ntheta_transform*nzeta_transform)
			B_dot_grad_zeta_pest_mnc(isurf,1) = sum(B_dot_grad_zeta_on_theta_pest_grid) &
				/ (ntheta_transform*nzeta_transform)
			Bsubu_pest_mnc(isurf,1) = sum(Bsubu_on_theta_pest_grid)/(ntheta_transform*nzeta_transform)
			Bsubv_pest_mnc(isurf,1) = sum(Bsubv_on_theta_pest_grid)/(ntheta_transform*nzeta_transform)
			factor = 2.0/(ntheta_transform*nzeta_transform)
			do imn=2,mnmax_transform
				angle = xm_transform(imn)*thetas_transform_2d - xn_transform(imn)*zetas_transform_2d
				sinangle = sin(angle)
				cosangle = cos(angle)
				factor2 = factor
				if (mod(ntheta_transform,2)==0 .and. xm_transform(imn)==(ntheta_transform/2)) then
					factor2 = factor2 / 2.0
				end if
				if (mod(nzeta_transform,2)==0 .and. xn_transform(imn)==(nfp*nzeta_transform/2)) then
					factor2 = factor2 / 2.0
				end if
				B_pest_mnc(isurf,imn) = sum(B_on_theta_pest_grid*cosangle)*factor2
				B_dot_grad_zeta_pest_mnc(isurf,imn) = sum(B_dot_grad_zeta_on_theta_pest_grid*cosangle)*factor2
				Bsubu_pest_mnc(isurf,imn) = sum(Bsubu_on_theta_pest_grid*cosangle)*factor2
				Bsubv_pest_mnc(isurf,imn) = sum(Bsubv_on_theta_pest_grid*cosangle)*factor2
			end do ! imn
			! Sanity tests - inverse FT
			do itheta=1,ntheta_transform
				do izeta=1,nzeta_transform
					B_on_theta_pest_grid_reconstructed(itheta,izeta) = &
						dot_product(B_pest_mnc(isurf,:),cos(xm_transform*thetas_transform(itheta) &
							- xn_transform*zetas_transform(izeta)))
					B_dot_grad_zeta_on_theta_pest_grid_reconstructed(itheta,izeta) = &
						dot_product(B_dot_grad_zeta_pest_mnc(isurf,:),cos(xm_transform*thetas_transform(itheta) &
							- xn_transform*zetas_transform(izeta)))
					Bsubu_on_theta_pest_grid_reconstructed(itheta,izeta) = &
						dot_product(Bsubu_pest_mnc(isurf,:),cos(xm_transform*thetas_transform(itheta) &
							- xn_transform*zetas_transform(izeta)))
					Bsubv_on_theta_pest_grid_reconstructed(itheta,izeta) = &
						dot_product(Bsubv_pest_mnc(isurf,:),cos(xm_transform*thetas_transform(itheta) &
							- xn_transform*zetas_transform(izeta)))
				end do
			end do
			! Compute error in FT
			if (any(abs(B_on_theta_pest_grid_reconstructed-B_on_theta_pest_grid)>1.0e-3)) then
				stop "Error in B_on_theta_pest_grid FT!"
			end if
			if (any(abs(B_dot_grad_zeta_on_theta_pest_grid_reconstructed-&
				B_dot_grad_zeta_on_theta_pest_grid)>1.0e-3)) then
				stop "Error in B_dot_grad_zeta_on_theta_pest_grid FT!"
			end if
			if (any(abs(Bsubu_on_theta_pest_grid_reconstructed-&
				Bsubu_on_theta_pest_grid)>1.0e-3)) then
				stop "Error in Bsubu_on_theta_pest_grid FT!"
			end if
			if (any(abs(Bsubv_on_theta_pest_grid_reconstructed-&
				Bsubv_on_theta_pest_grid)>1.0e-3)) then
				stop "Error in Bsubv_on_theta_pest_grid FT!"
			end if
		end do ! isurf

	end subroutine init_geometry_vmec

	function fzero_residual(theta_vmec_try)

		implicit none

		real(dp), intent(in) :: theta_vmec_try
		real(dp) :: fzero_residual
		integer :: imn

		fzero_residual = theta_vmec_try - theta_pest_target

		fzero_residual = fzero_residual &
			+ dot_product(lmns_slice,sin(xm_vmec*theta_vmec_try - xn_vmec*this_zeta))

	end function fzero_residual

end module geometry_mod
