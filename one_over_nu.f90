program one_over_nu

    use stel_kinds
    use splines_mod, only: init_spline
    use geometry_mod, only: init_geometry_boozer, init_geometry_vmec
    use grids_mod, only: init_grids
    use extrema_mod, only: find_extrema
    use input_mod, only: read_input, validate_input, geometry_option, verbose

    implicit none

    integer :: tic, toc, countrate
    real(dp) :: total_time

    call system_clock(tic,countrate)

    call read_input()
    call validate_input()
    call init_grids()
    if (geometry_option == 1) then
        call init_geometry_boozer()
    else
        call init_geometry_vmec()
    end if
    call find_extrema()
    call init_spline()
    call surf_loop()
    call system_clock(toc)
    total_time = real(toc-tic)/countrate
    if (verbose) then
        print *,"one_over_nu complete. Total time=",total_time,"sec."
    end if
    call write_output(total_time)

end program one_over_nu
