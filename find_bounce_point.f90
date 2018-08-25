! ===================================================
! Subroutine find_bounce_point
!
! This subroutine finds the bounce point using a
! Newton root solve. Bounce point must be bracketed
! by zeta_left and zeta_right. Much of this code
! has been taken from the rtsafe routine in
! Numerical Recipes. 
!
! Input:
! isurf				index in ssurf corresponding to calculation
!	zeta_left		left bound on zeta
! zeta_right	right bound on zeta
!	alpha0			field line label for calculation
! lambda			(scaled) lambda for calculation
!
! Output:
!	zeta_bounce_point		Computed bounce point
!
! ===================================================
subroutine find_bounce_point(isurf,zeta_left,zeta_right,alpha0,lambda,zeta_bounce_point)

	use stel_kinds
	use constants_mod
	use geometry_mod, only: iota
	use input_mod, only: niter_root, root_search_tolerance
	use splines_mod, only: compute_B_spline, compute_geometry_spline, compute_dBdzeta_spline

	implicit none

	real(dp), intent(in) :: zeta_left, zeta_right, alpha0, lambda
	real(dp), intent(out) :: zeta_bounce_point
	real(dp) :: modB_left, modB_right, fl, fh, xl, xh, dBdtheta, dBdzeta
	integer :: exit_code, j
	real(dp), dimension(Niter_root+1) :: zeta_search, f_search
	real(dp) :: dxold, dxrts, f, df, dx, modB, rts, temp
	real(dp), dimension(geometry_length) :: geometry
	integer :: isurf

	zeta_search = 0
	f_search = 0

	modB_left = compute_B_spline(isurf,alpha0+iota(isurf)*zeta_left,zeta_left)
	modB_right = compute_B_spline(isurf,alpha0+iota(isurf)*zeta_right,zeta_right)
	fl = 1-lambda*modB_left
	fh = 1-lambda*modB_right

	if ((fl > 0  .and. fh > 0) .or. (fl < 0 .and. fh < 0)) then
		print *,"Error! Root must be bracketed in newton solve."
		stop
	end if
	if (fl==0.0) then
		zeta_bounce_point = zeta_left
		return
	end if
	if (fh==0.0) then
		zeta_bounce_point = zeta_right
		return
	end if

	! Orient the search so that f(x1) < 0
	if (fl < 0) then
		xl = zeta_left
		xh = zeta_right
	else
		xh = zeta_left
		xl = zeta_right
	end if

	! Initialize guess for root, "stepsize before last", and last step
	rts = 0.5*(zeta_left+zeta_right)
	dxold = abs(zeta_left-zeta_right)
	dx = dxold
	! Evaluate function and derivative at initial guess
	geometry = compute_geometry_spline(isurf,alpha0+iota(isurf)*rts,rts)
	modB = geometry(B_index)
	f = 1-lambda*modB

	! Compute derivative of (1-lambda*B) wrt zeta
	dBdtheta = geometry(dBdtheta_index)
	dBdzeta = compute_dBdzeta_spline(isurf,alpha0+iota(isurf)*rts,rts)
	df = -lambda*(iota(isurf)*dBdtheta + dBdzeta)

	zeta_search(1) = rts
	f_search(1) = f

		exit_code = -1
    ! Loop over allowed iterations
    do j = 1, Niter_root
      ! Bisect if Newton out of range or not decreasing fast enough
      if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) .or. (abs(2*f) > abs(dxold*df))) then
        dxold = dx
        dx = 0.5*(xh-xl)
        rts = xl+dx
				! Change in root is negligible.
        if (xl == rts) then
          zeta_bounce_point = rts
          exit_code = 0
          exit
        end if
      ! Newton step acceptable. Take it.
      else
        dxold = dx
        dx = f/df
        temp = rts
        rts = rts - dx
        if (temp == rts) then
          zeta_bounce_point = rts
          exit_code = 0
          exit
        end if
      end if
      ! Convergence criterion
      if (abs(dx) < root_search_tolerance) then
        zeta_bounce_point = rts
        exit_code = 0
        exit
      end if
      ! The one new function evaluation per iteration
			geometry = compute_geometry_spline(isurf,alpha0+iota(isurf)*rts,rts)
      modB =  geometry(B_index)
			f = 1-lambda*modB

      zeta_search(j+1) = rts
      f_search(j+1) = f
			! Compute derivative of (1-lambda*B) wrt zeta
			dBdtheta = geometry(dBdtheta_index)
			dBdzeta = compute_dBdzeta_spline(isurf,alpha0+iota(isurf)*rts,rts)
      df = -lambda*(iota(isurf)*dBdtheta + dBdzeta)
      ! Maintain the bracket on the root
      if (f < 0) then
        xl = rts
      else
        xh = rts
      end if
    end do

	if (exit_code == -1) then
		print *,"*******************************************************************************"
		print *,"*******************************************************************************"
		print *,"The bounce point search did not converge within Niter_root iterations!"
		print *,"*******************************************************************************"
		print *,"*******************************************************************************"
		print *,"Here are the zetas we used: "
		print *,zeta_search
		print *,"Here are the values of 1-lambda*B: "
		print *,f_search
		stop
	end if

end subroutine find_bounce_point
