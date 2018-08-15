module geometry_mod

	use constants_mod
	use stel_kinds
	use input_mod, only: boozmn_filename, nsurf

	! Boozmn geometry items
	integer :: ns, nfp, nmodes
	real(dp), dimension(:), allocatable :: iota, Boozer_I, Boozer_G, xm, xn
	real(dp), dimension(:,:), allocatable :: bmnc
	! Value of s for nearest surface in boozmn_filename (used for calculation)
	real(dp), dimension(:), allocatable :: s_surf

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
!	compute_geometry	array of length geometry_length=4
!										This array can be indexed into
!										using B_index=1, dBdtheta=2,
!										dBdzeta=3 (defined in constants_mod).
!
! ===================================================
	function compute_geometry(isurf,theta,zeta)

		real(dp) :: theta,zeta
		real(dp), dimension(geometry_length) :: compute_geometry
		integer :: imn
		real(dp), dimension(nmodes) :: cosangle, sinangle
		integer, intent(in) :: isurf

		do imn=1,nmodes
			cosangle(imn) = cos(xm(imn)*theta-xn(imn)*zeta)
			sinangle(imn) = sin(xm(imn)*theta-xn(imn)*zeta)
		end do

		compute_geometry(B_index) = dot_product(bmnc(isurf,:),cosangle) ! B
		compute_geometry(dBdtheta_index) = -dot_product(bmnc(isurf,:),xm*sinangle) ! dBdtheta
		compute_geometry(dBdzeta_index) = dot_product(bmnc(isurf,:),xn*sinangle) ! dBdzeta

	end function

	real(dp) function compute_B(isurf,theta,zeta)
		real(dp) :: theta,zeta
		integer :: isurf

		compute_B = dot_product(bmnc(isurf,:),cos(xm*theta-xn*zeta))

	end function

! ===================================================
! subroutine init_geometry
!
!	This subroutine reads the boozmn input file
! to get the required geometric quantities on surfaces
! closest to the on requested surfaces in s_wish. This
! uses the read_boozer_mod from LIBSTELL.
!
! ===================================================
	subroutine init_geometry

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

	end subroutine init_geometry

end module geometry_mod
