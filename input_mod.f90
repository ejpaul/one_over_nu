module input_mod

	use stel_kinds
	use constants_mod

	! Input parameters for Newton minimization to find B extrema
	integer :: niter_newton = 10
	real(dp) :: tol_newton = 1e-3
	! Parameters for vmec coordinate transormation to straight field line theta
	real(dp) :: mpol_transform_refinement = 2.0
	real(dp) :: ntor_transform_refinement = 2.0
	real(dp) :: transform_abserr = 1.0e-10
	real(dp) :: transform_relerr = 1.0e-10
	! Grid dimensions
	integer :: nalpha = 100
	integer :: nlambda = 100
	! Geometry input
	! 1 = boozmn, 2 = VMEC
	integer :: geometry_option = 1
	integer :: nsurf ! length of s_wish
	character(len=200) :: boozmn_filename
	character(len=200) :: wout_filename
	! Input parameters for splines
	integer :: ntheta_spline = 200, nzeta_spline = 200
	! Resolution parameters
	integer :: nintegral = 10 ! Number of points of RK4 evaluation of bounce integral
	integer :: nwell = 100 ! Not really a resolution parameter - only controls output of classes
	! Stepsize for search for bounce points
	real(dp) :: Delta_zeta = 0.01

	! Maximum search in zeta for bounce point
	real(dp) :: max_search_in_zeta = 3000

	! Parameters for Newton root solve for bounce points
	real(dp) :: root_search_tolerance= 1e-3
	integer :: Niter_root = 10

	logical :: verbose = .true.

	! Option for computing particle flux for given profiles
	logical :: output_particle_flux = .false.
	! collision_species_option = 1 -> electron-ion collisions. q_e(2) used for charge of ion
	! collision_species_option = 2 -> self collisions. q_e(2) not used.
	integer :: collision_species_option = 1
	logical :: output_p_tensor = .false.
	logical :: output_J = .false.
	! Profiles of length s_wish_length
	! T_ev - Temperature (eletron volts)
	! n_m3 - density (m^{-3})
	! dlnTdpsi - temperature gradient (T^{-1} m^{-1})
	! dlnndpsi - density gradient (n^{-1} m^{-1})
	real(dp), dimension(s_wish_length) :: T_eV, dlnTdpsi, n_m3, dlnndpsi
	! m_kg - mass in kg
	! q_e - charge in elementary charge units (1) = primary (2) = secondary
	real(dp) :: m_kg
	integer, dimension(2) :: q_e

	namelist / one_over_nu_nml / nlambda, nalpha, nintegral, Delta_zeta, max_search_in_zeta, ntheta, nzeta, boozmn_filename, s_wish, nsurf, root_search_tolerance, Niter_root, nwell, niter_newton, tol_newton, verbose, ntheta_spline, nzeta_spline, output_particle_flux, T_ev, dTdpsi, n_m3, dlnndpsi, dlnTdpsi, m_kg, q_e, collision_species_option, output_p_tensor, geometry_option, ntor_transform_refinement, mpol_transform_refinement, transform_abserr, transform_relerr, wout_filename, output_J

	character(len=200) :: inputFilename

	! Desired surface for calculation
	real(dp), dimension(s_wish_length) :: s_wish

	contains

! ===================================================
! subroutine read_input
!
!	This subroutine reads the input namelist. Much of
! this code has been copied from regcoil_read_input.f90
!
! ===================================================
	subroutine read_input

		use constants_mod

		implicit none

		integer :: numargs, isurf
		integer :: fileUnit, didFileAccessWork, i
		integer, parameter :: uninitialized = -9999

		call getcarg(1, inputFilename, numargs)

		if (numargs<1) then
			 stop "One argument is required: the input namelist file, which must be named one_over_nu_in.XXXXX"
		end if
		if (numargs>1) then
			 print *,"WARNING: Arguments after the first will be ignored."
		end if
		if (inputFilename(1:15) .ne. "one_over_nu_in.") then
			 stop "Input file must be named one_over_nu_in.XXX for some extension XXX"
		end if

		fileUnit=11
		open(unit=fileUnit, file=inputFilename, action="read", status="old", iostat=didFileAccessWork)
		if (didFileAccessWork /= 0) then
			 print *,"Error opening input file ", trim(inputFilename)
			 stop
		else
			 read(fileUnit, nml=one_over_nu_nml, iostat=didFileAccessWork)
			 if (didFileAccessWork /= 0) then
					print *,"Error!  I was able to open the file ", trim(inputFilename), &
								 " but not read data from the one_over_nu_nml namelist in it."
					if (didFileAccessWork==-1) then
						 print *,"Make sure there is a carriage return after the / at the end of the namelist!"
					end if
					stop
			 end if
			 if (verbose) print *,"Successfully read parameters from one_over_nu_nml namelist in ", trim(inputFilename), "."
		end if
		close(unit = fileUnit)

		if (verbose) then
			 print *,"Input parameters:"
			 print "(a,i5)","nlambda =",nlambda
			 print "(a,i5)","nalpha  =",nalpha
			 print "(a,i5)","nintegral =",nintegral
			 print "(a,i5)","ntheta = ",ntheta
			 print "(a,i5)","nzeta = ",nzeta
			 print "(a,i5)","Niter_root = ",Niter_root
			 print "(a,i5)","niter_newton = ",niter_newton
			 print "(a,E15.7)","Delta_zeta = ", Delta_zeta
			 print "(a,E15.7)","root_search_tolerance", root_search_tolerance
			 print "(a,E15.7)","tol_newton = ", tol_newton
			 print "(a,i5)","nsurf = ",nsurf
			 print "(a,i5)","nwell = ",nwell
			 print "(a,i5)","ntheta_spline = ",ntheta_spline
			 print "(a,i5)","nzeta_spline = ",nzeta_spline
			 print "(a,i1)","geometry_option = ", geometry_option
			 print "(a,l)","output_J = ", output_J
			 if (geometry_option == 1) then
			 		print "(a,a)","boozmn_filename = ",boozmn_filename
			 else
					print "(a,a)","wout_filename = ",wout_filename
					print "(a,E15.7)","mpol_transform_refinement = ", mpol_transform_refinement
					print "(a,E15.7)","ntor_transform_refinement = ", ntor_transform_refinement
					print "(a,E15.7)","transform_abserr = ", transform_abserr
					print "(a,E15.7)","transform_relerr = ", transform_relerr
			 end if
			 print "(a,E15.7)","max_search_in_zeta = ", max_search_in_zeta
			 if (nsurf > 1) then
				 write(*,fmt = "(a,E15.7)", advance = "no") "s_wish = ", s_wish(1)
				 do isurf=2,nsurf-1
						write(*,fmt ="(E15.7)", advance = "no") s_wish(isurf)
				 end do
				 write(*,fmt ="(E15.7)", advance = "yes") s_wish(nsurf)
			 else
				 write(*,fmt = "(a,E15.7)", advance = "yes") "s_wish = ", s_wish(1)
			 end if
			 print "(a,l)","outut_particle_flux = ",output_particle_flux
			 if (output_particle_flux) then
					! Write T_ev
					if (nsurf > 1) then
						write(*,fmt = "(a,E15.7)", advance = "no") "T_ev = ", T_ev(1)
						do isurf=2,nsurf
							write(*,fmt ="(E15.7)", advance = "no") T_ev(isurf)
						end do
						write(*,fmt ="(E15.7)", advance = "yes") T_ev(nsurf)
					else
						write(*,fmt = "(a,E15.7)", advance = "yes") "T_ev = ", T_ev(1)
					end if
					! Write n_m3
					if (nsurf > 1) then
						write(*,fmt = "(a,E15.7)", advance = "no") "n_m3 = ", n_m3(1)
						do isurf=2,nsurf
							write(*,fmt ="(E15.7)", advance = "no") n_m3(isurf)
						end do
						write(*,fmt ="(E15.7)", advance = "yes") n_m3(nsurf)
					else
						write(*,fmt = "(a,E15.7)", advance = "yes") "n_m3 = ", n_m3(1)
					end if
					! Write dlnTdpsi
					if (nsurf > 1) then
						write(*,fmt = "(a,E15.7)", advance = "no") "dlnTdpsi = ", dlnTdpsi(1)
						do isurf=2,nsurf
							write(*,fmt ="(E15.7)", advance = "no") dlnTdpsi(isurf)
						end do
						write(*,fmt ="(E15.7)", advance = "yes") dlnTdpsi(nsurf)
					else
						write(*,fmt = "(a,E15.7)", advance = "yes") "dlnTdpsi = ", dlnTdpsi(1)
					end if
					! Write dlnndpsi
					if (nsurf > 1) then
						write(*,fmt = "(a,E15.7)", advance = "no") "dlnndpsi = ", dlnndpsi(1)
						do isurf=2,nsurf
							write(*,fmt ="(E15.7)", advance = "no") dlnndpsi(isurf)
						end do
						write(*,fmt ="(E15.7)", advance = "yes") dlnndpsi(nsurf)
					else
						write(*,fmt = "(a,E15.7)", advance = "yes") "dlnndpsi = ", dlnndpsi(1)
					end if
					print "(a,i1)","collision_species_option = ", collision_species_option
					print "(a,E15.7)","m_kg = ", m_kg
					if (collision_species_option==1) then
						print "(a,i1,i1)","q_e = ", q_e(1)
					else
						print "(a,i1,i1)","q_e = ", q_e(1),q_e(2)
					end if
			 end if
			 print "(a,l)","output_p_tensor = ", output_p_tensor
		end if

	end subroutine read_input

! ===================================================
! subroutine validate_input
!
!	This subroutine reads the checks the input parameters.
! Much of this has been copied from regcoil_validate_input.f90
!
! ===================================================
	subroutine validate_input

		use safe_open_mod

		implicit none

		integer :: isurf

		if (nlambda < 2) then
			stop "Error! nlambda must be >= 2."
		end if

		if (nalpha < 2) then
			stop "Error! nalpha must be >= 2."
		end if

		if (nintegral < 1) then
			stop "Error! nalpha must be >= 1."
		end if

		if (Delta_zeta <= 0) then
			stop "Error! Delta_zeta must be > 0."
		end if

		if (max_search_in_zeta <= 0) then
			stop "Error! max_search_in_zeta must be > 0."
		end if

		if (ntheta < 2) then
			stop "Error! ntheta must be >= 2."
		end if

		if (nzeta < 2) then
			stop "Error! nzeta must be >= 2."
		end if

		do isurf=1,nsurf
			if (s_wish(isurf)<0 .or. s_wish(isurf)>1) then
				stop "Error! s_wish must be >=0 and <=1."
			end if
		end do

		if (root_search_tolerance<=0 .or. root_search_tolerance>=1) then
			stop "Error! root_search_tolerance must be >0 and <1."
		end if

		if (Niter_root<1) then
			stop "Error! Niter_root must be >= 1."
		end if

		if (nwell<1) then
			stop "Error! nwell must be >= 1."
		end if

		if (niter_newton<1) then
			stop "Error! Niter_newton must be >= 1."
		end if

		if (tol_newton<=0 .or. tol_newton>=1) then
			stop "Error! tol_newton must be >0 and <1."
		end if

		if (ntheta_spline<=2) then
			stop "Error! ntheta_splien must be >=2."
		end if

		if (nzeta_spline<=2) then
			stop "Error! nzeta_spline must be >=2."
		end if

		if (collision_species_option > 2 .or. collision_species_option < 1) then
			stop "Erro! collision_species_option must be 0 or 1."
		end if

		if (output_particle_flux .and. collision_species_option==1) then
			if (verbose) then
				print *,"Electron-ion collision option chosen. q_e(2) will be used for secondary charge."
			end if
		end if
		if (output_particle_flux .and. collision_species_option==2) then
			if (verbose) then
				print *,"Self collision option chosen. q_e(2) will be ignored."
			end if
		end if
		if (geometry_option < 1 .or. geometry_option > 2) then
			stop "Error! Geometry_option must be 1 or 2."
		end if
		if (mpol_transform_refinement < 1 .or. ntor_transform_refinement < 1) then
			print *,"Warning. mpol_transform_refinement and mpol_transform_refinement should be > 1."
		end if
		if (transform_relerr < 0 .or. transform_relerr > 1) then
			print *,"transform_relerr should be between 0 and 1."
		end if
		if (transform_abserr < 0 .or. transform_abserr > 1) then
			print *,"transform_abserr should be between 0 and 1."
		end if

	end subroutine validate_input

end module input_mod

