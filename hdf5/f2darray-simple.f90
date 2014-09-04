program arraysimple
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
    rundata % filename = "data-simple-fort.h5"

    call get_options(rundata)

    allocate(dens(rundata%nx, rundata%ny))
    allocate(vel (2, rundata%nx, rundata%ny))
 
    call fillarray2d(dens, rundata % nx, rundata % ny)
    call fillarray3d(vel,  rundata % nx, rundata % ny)

    if ((rundata % nx)*(rundata % ny) < 200)  then
        call printarray2d(dens)
    endif
    call writehdf5file(rundata, dens, vel)

    deallocate(dens)
    deallocate(vel)


contains
    subroutine writehdf5file(rundata, dens, vel)
        use hdf5
        implicit none
        type(rundata_t), intent(IN) :: rundata
        double precision, intent(IN), dimension(:,:) :: dens
        double precision, intent(IN), dimension(:,:,:) :: vel

        integer(hid_t) :: file_id
        integer(hid_t) :: dens_space_id, vel_space_id
        integer(hid_t) :: dens_id, vel_id
        integer(hsize_t), dimension(2) :: densdims
        integer(hsize_t), dimension(3) :: veldims
        
        integer :: status

        ! first we have to open the FORTRAN interface.

        call h5open_f(status)
        ! create the file, check return code

        call h5fcreate_f(rundata%filename, H5F_ACC_TRUNC_F, file_id, status)
        if (status /= 0) then
            print *,'Could not open file ', rundata%filename
            return
        endif

        ! create the dataspaces corresponding to our variables
        densdims = (/ rundata % nx, rundata % ny /)
        call h5screate_simple_f(2, densdims, dens_space_id, status) 
        if (status /= 0) print *,'Could not create dens_space_id'

        veldims = (/ 2, rundata % nx, rundata % ny /)
        call h5screate_simple_f(3, veldims, vel_space_id, status) 
        if (status /= 0) print *,'Could not create vel_space_id'

        ! now that the dataspaces are defined, we can define variables on them

        call h5dcreate_f(file_id, "dens", H5T_IEEE_F64LE, dens_space_id, dens_id, status)
        if (status /= 0) print *,'Could not create dens_id'
        call h5dcreate_f(file_id, "vel" , H5T_IEEE_F64LE, vel_space_id,  vel_id, status)
        if (status /= 0) print *,'Could not create vel_id'

        ! Write the data.  We're writing it from memory, where it is saved 
        ! in NATIVE_DOUBLE format 
        call h5dwrite_f(dens_id, H5T_NATIVE_DOUBLE, dens, densdims, status)
        call h5dwrite_f(vel_id,  H5T_NATIVE_DOUBLE, vel,  veldims,  status)

        ! and now we're done
        call h5sclose_f(dens_space_id, status)
        call h5dclose_f(dens_id, status)
        call h5sclose_f(vel_space_id, status)
        call h5dclose_f(vel_id, status)
        call h5fclose_f(file_id, status)

        call h5close_f(status)
    end subroutine writehdf5file

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
            print *,'Usage: f2darray-simple [--help] [filename [nx [ny]]]'
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
end program arraysimple
