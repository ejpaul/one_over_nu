module input_mod

	use stel_kinds

	! Input parameters for Newton minimization to find B extrema
	integer :: niter_newton = 10
	real(dp) :: tol_newton = 1e-3
	! Grid dimensions
	integer :: nalpha = 100
	integer :: nlambda = 100
	! Geometry input
	integer :: nsurf ! length of s_wish
	character(len=200) :: boozmn_filename
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

	namelist / one_over_nu_nml / nlambda, nalpha, nintegral, Delta_zeta, max_search_in_zeta, ntheta, nzeta, boozmn_filename, s_wish, nsurf, root_search_tolerance, Niter_root, nwell, niter_newton, tol_newton, verbose, ntheta_spline, nzeta_spline

	character(len=200) :: inputFilename

	! Desired surface for calculation
	real(dp), dimension(:), allocatable :: s_wish

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

		allocate(s_wish(s_wish_length))

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
			 print "(a,i5)","   nlambda =",nlambda
			 print "(a,i5)","   nalpha  =",nalpha
			 print "(a,i5)","   nintegral =",nintegral
			 print "(a,i5)","   ntheta = ",ntheta
			 print "(a,i5)","   nzeta = ",nzeta
			 print "(a,i5)","   Niter_root = ",Niter_root
			 print "(a,i5)","   niter_newton = ",niter_newton
			 print "(a,E15.7)","   Delta_zeta = ", Delta_zeta
			 print "(a,E15.7)","   root_search_tolerance", root_search_tolerance
			 print "(a,E15.7)","   tol_newton = ", tol_newton
			 print "(a,i5)"," 	nsurf = ",nsurf
			 print "(a,i5)","   nwell = ",nwell
			 print "(a,i5)","   ntheta_spline = ",ntheta_spline
			 print "(a,i5)","   nzeta_spline = ",nzeta_spline
			 print "(a,a)","   boozmn_filename = ",boozmn_filename
			 print "(a,E15.7)","   max_search_in_zeta = ", max_search_in_zeta
			 write(*,fmt = "(a,E15.7)", advance = "no") "	 s_wish = ", s_wish(1)
			 do isurf=2,nsurf-1
					write(*,fmt ="(E15.7)", advance = "no") s_wish(isurf)
			 end do
			 if (nsurf>1) then
					write(*,fmt ="(E15.7)", advance = "yes") s_wish(nsurf)
			 end if
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

	end subroutine validate_input

end module input_mod

