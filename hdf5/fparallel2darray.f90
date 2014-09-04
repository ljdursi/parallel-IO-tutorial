program parallelarray
    use mpi
    implicit none

    type :: rundata_t
        integer :: globalnx, globalny
        integer :: localnx, localny
        integer :: npx, npy
        integer :: myx, myy
        integer :: rank, nprocs
        character(len=100) :: filename
    end type rundata_t

    double precision, allocatable :: dens(:,:)
    double precision, allocatable :: vel(:,:,:)
    type(rundata_t) ::  rundata

    integer :: ierr

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rundata % rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, rundata % nprocs, ierr)

    ! set default values

    call nearsquare(rundata % nprocs , rundata % npx, rundata % npy)
    rundata % globalnx = 100;
    rundata % globalny = 100;
    rundata % localnx = (rundata % globalnx / rundata % npx);
    rundata % localny = (rundata % globalny / rundata % npy);
    rundata % filename = "fparalleldata.h5"

    call get_options(rundata)

    allocate(dens(rundata % localnx, rundata % localny))
    allocate(vel (2, rundata % localnx, rundata % localny))
 
    call fillarray2d(rundata, dens)
    call fillarray3d(rundata, vel) 

    if ((rundata % localnx)*(rundata % localny) < 200)  then
        call printarray2d(rundata % rank,dens)
    endif
    call writehdf5file(rundata, dens, vel)

    call MPI_Finalize(ierr)

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
        integer(hid_t) :: mem_dens_space_id, mem_vel_space_id
        integer(hid_t) :: dens_id, vel_id
        integer(hid_t) :: arr_group_id, attr_space_id, attr_id
        integer(hid_t) :: other_group_id
        integer(hid_t) :: timestep_id, timestep_space_id
        integer(hid_t) :: comptime_id, comptime_space_id
        integer(hid_t) :: author_id, author_space_id, author_type_id
        integer(hsize_t), dimension(2) :: densdims
        integer(hsize_t), dimension(3) :: veldims
        integer(hsize_t), dimension(2) :: locdensdims
        integer(hsize_t), dimension(3) :: locveldims
        integer(hsize_t), dimension(1) :: dummydims
        integer(hsize_t)  :: strlen
        integer, dimension(3) :: ddmmyy
        integer :: yyyymm
        
        integer :: status

        character(len=100):: authorname
        integer           :: timestep
        real              :: comptime

        integer :: info
        integer(hid_t) :: fap_id, dist_id
        integer(hsize_t), dimension(2) :: offsets, counts, strides, blocks
        integer(hsize_t), dimension(3) :: voffsets, vcounts, vstrides, vblocks
    
        dummydims = (/1/)


        ! first we have to open the FORTRAN interface.

        call h5open_f(status)

        ! An example of setting the MPI-IO hints for better performance on GPFS 
        call MPI_Info_create(info, status)
        call MPI_Info_set(info,"IBM_largeblock_io","true", status)

        ! A file accesss parameter:
        call H5Pcreate_f(H5P_FILE_ACCESS_F, fap_id, status)

        ! Include the file access property with IBM hint 
        call H5Pset_fapl_mpio_f(fap_id, MPI_COMM_WORLD, info, status)

        ! Create a new file - truncate anything existing, use these new properties
        call H5Fcreate_f(rundata%filename, H5F_ACC_TRUNC_F, file_id, status, access_prp=fap_id)
        if (status /= 0) then
            print *,'Could not open file ', rundata%filename
            return
        endif

        ! create a new group within this file
        call h5gcreate_f(file_id, '/ArrayData', arr_group_id, status)

        ! Give this group an attribute listing the time of calculation 
        call idate(ddmmyy)
        yyyymm = ddmmyy(3)*100+ddmmyy(2)

        call h5screate_f(H5S_SCALAR_F, attr_space_id, status)
        call h5acreate_f(arr_group_id, 'Calculated on (YYYYMM)', H5T_STD_U32LE, attr_space_id, attr_id, status)
        call h5awrite_f (attr_id, H5T_NATIVE_INTEGER, yyyymm, dummydims, status)
 
        ! we're done with the attribute and its space now
        call h5aclose_f(attr_id, status)
        call h5sclose_f(attr_space_id, status)

        ! create the dataspaces corresponding to our variables
        ! note that the file space contains the *global* view
        densdims = (/ rundata % globalnx, rundata % globalny /)
        call h5screate_simple_f(2, densdims, dens_space_id, status) 
        if (status /= 0) print *,'Could not create dens_space_id'

        veldims = (/ 2, rundata % globalnx, rundata % globalny /)
        call h5screate_simple_f(3, veldims, vel_space_id, status) 
        if (status /= 0) print *,'Could not create vel_space_id'


        ! now that the dataspaces are defined, we can define variables on them

        ! Note that we can create the file explicitly from the root group within the file:
        call h5dcreate_f(file_id, "/ArrayData/dens", H5T_IEEE_F64LE, dens_space_id, dens_id, status)
        if (status /= 0) print *,'Could not create dens_id'

        ! or simply within the group_id of the group we want it in...
        call h5dcreate_f(arr_group_id, "vel" , H5T_IEEE_F64LE, vel_space_id,  vel_id, status)
        if (status /= 0) print *,'Could not create vel_id'

        ! A transfer descriptor, describing how transfer is done:
        call H5Pcreate_f(H5P_DATASET_XFER_F, dist_id, status)

        ! we'll be transferring data with mpiio, collectively 
        call H5Pset_dxpl_mpio_f(dist_id, H5FD_MPIO_COLLECTIVE_F, status)


        !
        ! Now we have to figure out the `hyperslab' within the global
        ! data that corresponds to our local data.
        !
        ! Hyperslabs are described by an array of counts, strides, offsets,
        ! and block sizes.
        !
        !       |-offx--|
        !       +-------|----|-------+   -+-
        !       |                    |    |
        !       |                    |   offy
        !       |                    |    |
        !       -       +----+       -   -+-
        !       |       |    |       |    |
        !       |       |    |       |  localny
        !       |       |    |       |    |
        !       -       +----+       -   -+-
        !       |                    |
        !       |                    |
        !       +-------|----|-------+
        !               localnx
        !
        !  In this case the blocksizes are (localnx,localny) and the offsets are 
        !  (offx,offy) = ((myx)/nxp*globalnx, (myy/nyp)*globalny)
        !

        offsets(1) = (rundata % globalnx/rundata % npx)*rundata % myx 
        offsets(2) = (rundata % globalny/rundata % npy)*rundata % myy
        
        blocks(1)  = rundata % localnx
        blocks(2)  = rundata % localny

        counts(1)  = 1
        counts(2)  = 1

        strides(1)  = 1
        strides(2)  = 1

        ! select this subset of the density variable's space in the file 
        call H5Sselect_hyperslab_f(dens_space_id,H5S_SELECT_SET_F, offsets, &
                                   counts, status, strides, blocks)

        ! For the velocities, it's the same thing but there's a count of two,
        ! (one for each velocity component) */

        voffsets(1) = 0
        voffsets(2) = (rundata % globalnx/rundata % npx)*rundata % myx
        voffsets(3) = (rundata % globalny/rundata % npy)*rundata % myy
        
        vblocks(1)  = 1
        vblocks(2)  = rundata % localnx
        vblocks(3)  = rundata % localny

        vstrides(1)  = 1
        vstrides (2) = 1
        vstrides (3) = 1

        vcounts (1) = 2
        vcounts (2) = 1
        vcounts (3) = 1

        ! select this subset of the density variable's space in the file 
        call H5Sselect_hyperslab_f(vel_space_id,H5S_SELECT_SET_F, voffsets, &
                                   vcounts, status, vstrides, vblocks)


        ! Write the data.  We're writing it from memory, where it is saved 
        ! in NATIVE_DOUBLE format 
        locdensdims = (/ rundata % localnx, rundata % localny /)
        call H5Screate_simple_f(2, locdensdims, mem_dens_space_id, status)
        call h5dwrite_f(dens_id, H5T_NATIVE_DOUBLE, dens, locdensdims, status, &
                        mem_space_id=mem_dens_space_id, file_space_id=dens_space_id, &
                        xfer_prp=dist_id)

        locveldims = (/ 2, rundata % localnx, rundata % localny /)
        call H5Screate_simple_f(3, locveldims,  mem_vel_space_id, status)
        call h5dwrite_f(vel_id,  H5T_NATIVE_DOUBLE, vel,  locveldims,  status, &
                        mem_space_id=mem_vel_space_id, file_space_id=vel_space_id, &
                        xfer_prp=dist_id)

        ! We'll create another group for related info and put some things in there 

        authorname="Jonathan Dursi"
        timestep=13
        comptime=81.773

        ! create group 
        call h5gcreate_f(file_id, '/OtherStuff', other_group_id, status)

        ! scalar space, data for integer timestep 
        call h5screate_f(H5S_SCALAR_F, timestep_space_id, status)
        call h5dcreate_f(other_group_id, 'Timestep', H5T_STD_U32LE, timestep_space_id, timestep_id, status)
        call h5dwrite_f(timestep_id, H5T_NATIVE_INTEGER, timestep, dummydims, status)
        call h5dclose_f(timestep_id, status)
        call h5sclose_f(timestep_space_id, status)
       
        ! scalar space, data for floating compute time 
        call h5screate_f(H5S_SCALAR_F, comptime_space_id, status)
        call h5dcreate_f(other_group_id, 'Compute Time', H5T_IEEE_F32LE, comptime_space_id, comptime_id, status)
        call h5dwrite_f(comptime_id, H5T_NATIVE_REAL, comptime, dummydims, status)
        call h5dclose_f(comptime_id, status)
        call h5sclose_f(comptime_space_id, status)
       
        ! scalar space, data for string author name
        call h5screate_f(H5S_SCALAR_F, author_space_id, status)
        ! copy character type..
        call h5tcopy_f(H5T_NATIVE_CHARACTER, author_type_id, status)   
        ! and make it longer
        strlen = len_trim(authorname)
        call h5tset_size_f(author_type_id, strlen, status)
        call h5dcreate_f(other_group_id, 'Simulator Name', author_type_id, author_space_id, author_id, status)
        call h5dwrite_f(author_id, author_type_id, trim(authorname), dummydims, status)
        call h5dclose_f(author_id, status)
        call h5sclose_f(author_space_id, status)
        call h5tclose_f(author_type_id, status)
       
        call h5gclose_f(other_group_id, status)

        ! and now we're done

        call h5pclose_f(fap_id, status)
        call h5pclose_f(dist_id, status)
        call MPI_Info_free(info, status)

        call h5sclose_f(dens_space_id, status)
        call h5dclose_f(dens_id, status)
        call h5sclose_f(vel_space_id, status)
        call h5dclose_f(vel_id, status)
        call h5gclose_f(arr_group_id, status)
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

        if (nargs /= 0) then
            call getarg(1,arg)
            if ((trim(arg) == "-h") .or. (trim(arg) == "--help")) then
                print *,'Usage: parallelf2darray [--help] [filename [npx [npy [nx [ny]]]]]'
                print *,'       where filename is output filename,' 
                print *,'       npx, npy are number of prcessors along x, y directions, and'
                print *,'       nx, ny are total number of points in x and y directions.'
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
        endif
        

        if (nargs >= 2) then
            call getarg(2,arg)
            read(arg,*) iarg
            if ((iarg < 1) .or. (iarg > rundata % nprocs)) then
                print *,'Cannot use number of processors in X dimension ', arg, '; skipping.'
                print *,'Using ', rundata % npx, ' instead. '
            else  if (modulo(rundata % nprocs, iarg) /= 0) then
                print *,'Number of X processors ', iarg, ' does not divide ', rundata % nprocs,';  skipping.'
                print *,'Using ', rundata % npx, ' instead. '
            else
                rundata % npx = iarg
                rundata % npy = rundata % nprocs / iarg
            endif
        endif

        if (nargs >= 3) then
            call getarg(3,arg)
            read(arg,*) iarg
            if ((iarg < 1) .or. (iarg > rundata % nprocs)) then
                print *,'Cannot use number of processors in Y dimension ', arg, '; skipping.'
                print *,'Using ', rundata % npy, ' instead. '
            else  if (modulo(rundata % nprocs, iarg) /= 0) then
                print *,'Number of Y processors ', iarg, ' does not divide ', rundata % nprocs,' ; skipping.'
                print *,'Using ', rundata % npy, ' instead. '
            else
                rundata % npy = iarg
                rundata % npx = rundata % nprocs / iarg
            endif
        endif

        if (nargs >= 4) then
            call getarg(4,arg)
            read(arg,*) iarg
            if ((iarg < rundata % npx) .or. (iarg > maxpts)) then
                print *,'Cannot use number of x-points ', arg, '; skipping.'
            else
                rundata % globalnx = iarg
            endif
        endif

        if (nargs >= 5) then
            call getarg(5,arg)
            read(arg,*) iarg
            if ((iarg < rundata % npy) .or. (iarg > maxpts)) then
                print *,'Cannot use number of y-points ', arg, '; skipping.'
            else
                rundata % globalny = iarg
            endif
        endif

        ! figure out where we are in the 2d grid of processors 

        rundata % myy = rundata % rank / (rundata % npx)
        rundata % myx = modulo(rundata % rank , rundata % npx)

        ! last row/column gets any extra or fewer points to make things
        ! work out:

        rundata % localnx = rundata % globalnx / (rundata % npx)
        rundata % localny = rundata % globalny / (rundata % npy)

        if (rundata % myx == (rundata % npx - 1)) then
            rundata % localnx = rundata % globalnx - (rundata%npx - 1)*(rundata%localnx)
        endif
        if (rundata % myy == (rundata % npy - 1)) then
            rundata % localny = rundata % globalny - (rundata%npy - 1)*(rundata%localny)
        endif
 
        print '(A,I3,A,I2,A,I2,A,I3,A,I3,A,I3,A,I3,A)', &
                     '[',rundata % rank, '] gets (', rundata % myx, ', ', &
                     rundata % myy, '): local points = (', rundata % localnx, &
                     ',', rundata % localny, '); global points = (', &
                     rundata % globalnx, ',', rundata % globalny, ').'
 
    end subroutine get_options

    subroutine fillarray2d(rundata, dens)
        implicit none
        type(rundata_t), intent(in) :: rundata
        double precision, intent(out), dimension(:,:) :: dens
        
        integer :: i,j
        double precision :: sigma
        double precision, dimension(rundata % localnx) :: rx2, r2
    
        integer :: gnx, gny, nx, ny, npx, npy, myx, myy
        integer :: startx, starty

        gnx = rundata % globalnx
        gny = rundata % globalny
        nx  = rundata % localnx
        ny  = rundata % localny
        npx = rundata % npx
        npy = rundata % npy
        myx = rundata % myx
        myy = rundata % myy

        startx = (gnx/npx)*myx
        starty = (gny/npy)*myy

        sigma = gnx/4.

        rx2 = (/ (((i+startx)-gnx/2.)*((i+startx)-gnx/2.), i=1,nx) /)
        do j=1,ny
            r2 = rx2 + ((j+starty)-gny/2.)*((j+starty)-gny/2.)
            dens(:,j) = 1. + 4.*exp(-r2/(2.*sigma*sigma))
        enddo
    end subroutine fillarray2d

    subroutine fillarray3d(rundata, vel)
        implicit none
        type(rundata_t), intent(in) :: rundata
        double precision, intent(out), dimension(:,:,:) :: vel
        
        integer :: i,j
        double precision :: sigma
        double precision, dimension(rundata % localnx) :: rx2, r2

        integer :: gnx, gny, nx, ny, npx, npy, myx, myy
        integer :: startx, starty

        gnx = rundata % globalnx
        gny = rundata % globalny
        nx  = rundata % localnx
        ny  = rundata % localny
        npx = rundata % npx
        npy = rundata % npy
        myx = rundata % myx
        myy = rundata % myy

        startx = (gnx/npx)*myx
        starty = (gny/npy)*myy

        sigma = nx/4.

        rx2 = (/ (((i+startx)-gnx/2.)*((i+startx)-gnx/2.), i=1,nx) /)
        do j=1,ny
            r2 = rx2 + ((j+starty)-gny/2.)*((j+starty)-gny/2.)
            do i=1,nx 
                vel(1,i,j) = 4.*exp(-r2(i)/(2.*sigma*sigma))*((j+starty)-gny/2.)
                vel(2,i,j) =-4.*exp(-r2(i)/(2.*sigma*sigma))*((i+startx)-gnx/2.)
            enddo
        enddo
    end subroutine fillarray3d

    subroutine printarray2d(rank, dens)
        implicit none
        integer, intent(in) :: rank
        double precision, intent(in), dimension(:,:) :: dens

        print *, '[',rank,']', dens
    end subroutine printarray2d

    subroutine nearsquare(nprocs, npx, npy)
        implicit none
        integer, intent(IN)  :: nprocs
        integer, intent(OUT) :: npx, npy

        integer :: sq, n, m
        logical :: first

        sq = ceiling(sqrt(nprocs*1.d0))
        first = .true.
        if (sq*sq == nprocs) then
            npx = sq
            npy = sq
        else
            do n=sq+1,1,-1
               if (first .and. (modulo(nprocs,n) == 0)) then
                   npx = n
                   npy = nprocs/n
                   first = .false.
               endif
            enddo
        endif
    end subroutine nearsquare
end program parallelarray
