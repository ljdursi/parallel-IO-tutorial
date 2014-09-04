program readarraysimple
    implicit none

    type :: rundata_t
        integer :: nx, ny, nvelcomp
        character(len=100) :: filename
    end type rundata_t

    double precision, allocatable :: dens(:,:)
    double precision, allocatable :: vel(:,:,:)
    type(rundata_t) ::  rundata

    ! set default values

    rundata % filename = "data-simple-fort.h5"

    call get_options(rundata)
    call readhdf5file(rundata, dens, vel)

    print *, 'Read in: Dens(',rundata%nx,',',rundata%ny,')'
    print *, '         Vel (',rundata%nvelcomp,',',rundata%nx,',',rundata%ny,')'
    if ((rundata % nx)*(rundata % ny) < 200)  then
        call printarray2d(dens)
    endif

    deallocate(dens)
    deallocate(vel)


contains
    subroutine readhdf5file(rundata, dens, vel)
        use hdf5
        implicit none
        type(rundata_t), intent(INOUT) :: rundata
        double precision, allocatable, intent(OUT), dimension(:,:) :: dens
        double precision, allocatable, intent(OUT), dimension(:,:,:) :: vel

        integer(hid_t) :: file_id
        integer(hid_t) :: dens_space_id, vel_space_id
        integer(hid_t) :: dens_id, vel_id
        integer(hid_t) :: dtype
        integer(hsize_t), dimension(2) :: densdims
        integer(hsize_t), dimension(3) :: veldims
        integer(hsize_t), dimension(2) :: maxdensdims
        integer(hsize_t), dimension(3) :: maxveldims
        
        integer :: status
        integer :: ndims

        ! first we have to open the FORTRAN interface.
        call h5open_f(status)

        ! open the file, check return code
        call h5fopen_f(rundata%filename, H5F_ACC_RDONLY_F, file_id, status)
        if (status /= 0) then
            print *,'Could not open file ', rundata%filename
            return
        endif

        ! open the datasets
        call h5dopen_f(file_id, 'dens', dens_id, status)
        if (status /= 0) print *,'Could not open dens_id'
        call h5dopen_f(file_id, 'vel' , vel_id, status)
        if (status /= 0) print *,'Could not open vel_id'

        ! get the dataspaces each variable inhabit
        call h5dget_space_f(dens_id, dens_space_id, status) 
        if (status /= 0) print *,'Could not get dens_space_id'
        call h5dget_space_f(vel_id, vel_space_id, status) 
        if (status /= 0) print *,'Could not get vel_space_id'

        call h5sget_simple_extent_dims_f(dens_space_id, densdims, maxdensdims, ndims)
        if (ndims /= 2) print *,'Unexpected number of dens dimensions'
        call h5sget_simple_extent_dims_f(vel_space_id,  veldims,  maxveldims, ndims)
        if (ndims /= 0) print *,'Unexpected number of vel dimensions'

        ! now we can allocate the arrays
        allocate(dens(densdims(1), densdims(2)))
        allocate(vel(veldims(1), veldims(2), veldims(3)))

        rundata % nx = densdims(1)
        rundata % ny = densdims(2)
        rundata % nvelcomp = veldims(1)
    
        ! Read the data.  We're writing it from memory, where it is saved 
        ! in NATIVE_DOUBLE format 
        call h5dread_f(dens_id, H5T_NATIVE_DOUBLE, dens, densdims, status)
        call h5dread_f(vel_id,  H5T_NATIVE_DOUBLE, vel,  veldims,  status)

        ! and now we're done
        call h5sclose_f(dens_space_id, status)
        call h5dclose_f(dens_id, status)
        call h5sclose_f(vel_space_id, status)
        call h5dclose_f(vel_id, status)
        call h5fclose_f(file_id, status)

        call h5close_f(status)
    end subroutine readhdf5file

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
end program readarraysimple
