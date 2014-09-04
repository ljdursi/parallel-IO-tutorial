program array
    implicit none

    type :: rundata_t
        integer :: nx, ny
        character(len=100) :: filename
    end type rundata_t

    double precision, allocatable :: dens(:,:)
    double precision, allocatable :: vel(:,:,:)
    type(rundata_t) ::  rundata

    ! set default values

    rundata % nx = 100
    rundata % ny = 100
    rundata % filename = "data-fort.nc"

    call get_options(rundata)

    allocate(dens(rundata%nx, rundata%ny))
    allocate(vel (2, rundata%nx, rundata%ny))
 
    call fillarray2d(dens, rundata % nx, rundata % ny)
    call fillarray3d(vel,  rundata % nx, rundata % ny)

    if ((rundata % nx)*(rundata % ny) < 200)  then
        call printarray2d(dens)
    endif
    call writenetcdffile(rundata, dens, vel)

    deallocate(dens)
    deallocate(vel)


contains
    subroutine writenetcdffile(rundata, dens, vel)
        use netcdf
        implicit none
        type(rundata_t), intent(IN) :: rundata
        double precision, intent(IN), dimension(:,:) :: dens
        double precision, intent(IN), dimension(:,:,:) :: vel

        integer :: file_id, xdim_id, ydim_id, vcomp_id
        integer :: xcoord_id, ycoord_id
        integer :: dens_id, vel_id
        integer, dimension(2) :: densdims
        integer, dimension(3) :: veldims
        real, allocatable, dimension(:) :: x, y
        character(len=*), parameter :: coordunit = 'cm'
        character(len=*), parameter :: densunit  = 'g/cm^3'
        character(len=*), parameter :: velunit   = 'cm/s'
        
        integer :: i
        integer :: status

        i = 1
        ! create the file, check return code

        status = nf90_create(path=rundata%filename, cmode=NF90_CLOBBER, ncid=file_id)
        if (status /= NF90_NOERR) then
            print *,'Could not open file ', rundata%filename
            return
        endif

        ! define the dimensions

        status = nf90_def_dim(file_id, 'X', rundata%nx, xdim_id)
        status = nf90_def_dim(file_id, 'Y', rundata%ny, ydim_id)
        status = nf90_def_dim(file_id, 'velocity components', 2, vcomp_id)

        ! now that the dimensions are defined, we can define variables on them,...

        status = nf90_def_var(file_id, 'X coordinate', NF90_REAL, (/ xdim_id /), xcoord_id)
        status = nf90_def_var(file_id, 'Y coordinate', NF90_REAL, (/ ydim_id /), ycoord_id)

        densdims = (/ xdim_id, ydim_id /)
        veldims = (/ vcomp_id, xdim_id, ydim_id /)

        status = nf90_def_var(file_id, 'Density',  NF90_DOUBLE, densdims, dens_id)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Dens'
        status = nf90_def_var(file_id, 'Velocity', NF90_DOUBLE, veldims,  vel_id)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Vel'

         ! ...and assign units to them as an attribute 
         status = nf90_put_att(file_id, xcoord_id, "units", coordunit)
         status = nf90_put_att(file_id, ycoord_id, "units", coordunit)
         status = nf90_put_att(file_id, dens_id,   "units", densunit)
         status = nf90_put_att(file_id, vel_id,    "units", velunit)

        ! done defining

        status = nf90_enddef(file_id)

        ! Write out the values
        allocate(x(rundata%nx), y(rundata%ny))
        x = (/ (i-rundata%nx/2., i=1,rundata%nx) /)
        y = (/ (i-rundata%ny/2., i=1,rundata%ny) /)
        
        status = nf90_put_var(file_id, xcoord_id, x)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' X coord'
        status = nf90_put_var(file_id, ycoord_id, y)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Y coord'
        status = nf90_put_var(file_id, dens_id, dens)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Dens'
        status = nf90_put_var(file_id, vel_id,  vel)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Vel'

        status = nf90_close(file_id)
    
        deallocate(x,y)
    return
    end subroutine writenetcdffile

    subroutine get_options(rundata)
        implicit none
        type(rundata_t), intent(inout) :: rundata
        integer :: nargs
        character(100) :: arg
        integer, parameter :: maxpts = 500
        integer :: err
        integer :: iarg

        nargs = iargc()
        if (nargs == 0) return

        call getarg(1,arg)
        if ((trim(arg) == "-h") .or. (trim(arg) == "--help")) then
            print *,'Usage: f2darray [--help] [filename [nx [ny]]]'
            print *,'       where filename is output filename, and'
            print *,'       nx, ny are number of points in x and y directions.'
            call exit(0)
        endif
    
        !! at least one var - must be filename
        open(unit=14,file=arg,status='new',iostat=err)
        if (err /= 0) then
            print *,'Could not open file ', arg
            print *,'Exiting'
            call exit(1)
        endif
        close(unit=14)
        
        rundata % filename = arg

        if (nargs >= 2) then
            call getarg(2,arg)
            read(arg,*) iarg
            if ((iarg < 5) .or. (iarg > maxpts)) then
                print *,'Cannot use number of x-points ', arg, '; skipping.'
            else
                rundata % nx = iarg
            endif
        endif

        if (nargs >= 3) then
            call getarg(3,arg)
            read(arg,*) iarg
            if ((iarg < 5) .or. (iarg > maxpts)) then
                print *,'Cannot use number of y-points ', arg, '; skipping.'
            else
                rundata % ny = iarg
            endif
        endif

    end subroutine get_options

    subroutine fillarray2d(dens, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny 
        double precision, intent(out), dimension(:,:) :: dens
        
        integer :: i,j
        double precision :: sigma
        double precision, dimension(nx) :: rx2, r2

        sigma = nx/4.
        i = 1

        rx2 = (/ ((i-nx/2.)*(i-nx/2.), i=1,nx) /)
        do j=1,ny
            r2 = rx2 + (j-ny/2.)*(j-ny/2.)
            dens(:,j) = 1. + 4.*exp(-r2/(2.*sigma*sigma))
        enddo
    end subroutine fillarray2d

    subroutine fillarray3d(vel, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny 
        double precision, intent(out), dimension(:,:,:) :: vel
        
        integer :: i,j
        double precision :: sigma
        double precision, dimension(nx) :: rx2, r2

        sigma = nx/4.
        i = 1

        rx2 = (/ ((i-nx/2.)*(i-nx/2.), i=1,nx) /)
        do j=1,ny
            r2 = rx2 + (j-ny/2.)*(j-ny/2.)
            do i=1,nx 
                vel(1,i,j) = 4.*exp(-r2(i)/(2.*sigma*sigma))*(j-ny/2.)
                vel(2,i,j) =-4.*exp(-r2(i)/(2.*sigma*sigma))*(i-nx/2.)
            enddo
        enddo
    end subroutine fillarray3d

    subroutine printarray2d(dens)
        implicit none
        double precision, intent(in), dimension(:,:) :: dens

        print *, dens
    end subroutine printarray2d
end program array
