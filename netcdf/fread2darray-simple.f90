program readarraysimple
    implicit none

    type :: rundata_t
        integer :: nx, ny
        integer :: nvelcomp
        character(len=100) :: filename
    end type rundata_t

    double precision, allocatable :: dens(:,:)
    double precision, allocatable :: vel(:,:,:)
    type(rundata_t) ::  rundata

    ! set default values

    rundata % nx = 100
    rundata % ny = 100
    rundata % filename = "data-simple-fort.nc"

    call get_options(rundata)
    call readnetcdffile(rundata, dens, vel)

    print *, 'Read file ', rundata % filename
    print *, ' Density [', rundata % nx, '][ ', rundata % ny, '] '
    print *, ' Velocity[', rundata % nvelcomp, '][', rundata % nx, '][ ', rundata % ny, '] '

    if ((rundata % nx)*(rundata % ny) < 200)  then
        call printarray2d(dens)
    endif

    deallocate(dens)
    deallocate(vel)


contains
    subroutine readnetcdffile(rundata, dens, vel)
        use netcdf
        implicit none
        type(rundata_t), intent(INOUT) :: rundata
        double precision, allocatable, dimension(:,:), intent(OUT) :: dens
        double precision, allocatable, dimension(:,:,:), intent(OUT) :: vel

        integer :: file_id, xdim_id, ydim_id, vcomp_id
        integer :: dens_id, vel_id
        
        integer :: status

        ! create the file, check return code

        status = nf90_open(path=rundata%filename, mode=NF90_NOWRITE, ncid=file_id)
        if (status /= NF90_NOERR) then
            print *,'Could not open file ', rundata%filename
            return
        endif

        ! find the dimensions

        status = nf90_inq_dimid(file_id, 'X', xdim_id)
        status = nf90_inq_dimid(file_id, 'Y', ydim_id)
        status = nf90_inq_dimid(file_id, 'velocity components', vcomp_id)

        ! find the dimension lengths

        status = nf90_inquire_dimension(file_id, xdim_id, len = rundata % nx)
        status = nf90_inquire_dimension(file_id, ydim_id, len = rundata % ny )
        status = nf90_inquire_dimension(file_id, vcomp_id, len = rundata % nvelcomp)
        
        ! now we can allocate variable sizes

        allocate(dens(rundata%nx, rundata%ny)) 
        allocate(vel (rundata%nvelcomp, rundata%nx, rundata%ny)) 

        status = nf90_inq_varid(file_id, 'Density',  dens_id)
        status = nf90_inq_varid(file_id, 'Velocity', vel_id)

        ! read in the values
        status = nf90_get_var(file_id, dens_id, dens)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Dens'
        status = nf90_get_var(file_id, vel_id,  vel)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Vel'

        status = nf90_close(file_id)
    end subroutine readnetcdffile

    subroutine get_options(rundata)
        implicit none
        type(rundata_t), intent(inout) :: rundata
        integer :: nargs
        character(100) :: arg
        integer, parameter :: maxpts = 500
        integer :: err

        nargs = iargc()
        if (nargs == 0) return

        call getarg(1,arg)
        if ((trim(arg) == "-h") .or. (trim(arg) == "--help")) then
            print *,'Usage: fread2darray-simple [--help] [filename]'
            print *,'       where filename is input file name.'
            call exit(0)
        endif
    
        !! at least one var - must be filename
        open(unit=14,file=arg,status='old',iostat=err)
        if (err /= 0) then
            print *,'Could not open file ', arg
            print *,'Exiting'
            call exit(1)
        endif
        close(unit=14)
        
        rundata % filename = arg

    end subroutine get_options

    subroutine printarray2d(dens)
        implicit none
        double precision, intent(in), dimension(:,:) :: dens

        print *, dens
    end subroutine printarray2d
end program readarraysimple
