program one_over_nu

	use stel_kinds
	use splines_mod, only: init_spline
	use geometry_mod, only: init_geometry
	use grids_mod, only: init_grids
	use extrema_mod, only: find_extrema
	use input_mod, only: read_input, validate_input

  implicit none

  integer :: tic, toc, countrate
	real(dp) :: total_time

  call system_clock(tic,countrate)

  call read_input()
  call validate_input()
	call init_grids()
  call init_geometry()
	call find_extrema()
  call init_spline()
	call surf_loop()

	call system_clock(toc)
	total_time = real(toc-tic)/countrate
	print *,"one_over_nu complete. Total time=",total_time,"sec."

	call write_output(total_time)

end program one_over_nu
