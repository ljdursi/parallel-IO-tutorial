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
    rundata % filename = "fparalleldata.nc"

    call get_options(rundata)

    allocate(dens(rundata % localnx, rundata % localny))
    allocate(vel (2, rundata % localnx, rundata % localny))
 
    call fillarray2d(rundata, dens)
    call fillarray3d(rundata, vel) 

    if ((rundata % localnx)*(rundata % localny) < 200)  then
        call printarray2d(rundata % rank,dens)
    endif
    call writenetcdffile(rundata, dens, vel)

    call MPI_Finalize(ierr)

    deallocate(dens)
    deallocate(vel)


contains
    subroutine writenetcdffile(rundata, dens, vel)
        use mpi
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
        integer :: info
        integer :: mode_flag
    
        integer, dimension(2) :: densstarts
        integer, dimension(2) :: denscounts
        integer, dimension(3) :: velstarts
        integer, dimension(3) :: velcounts

        ! create the file, check return code

        call MPI_Info_create(info, status)
        call MPI_Info_set(info,"IBM_largeblock_io","true", status)
   
        mode_flag = IOR(NF90_MPIIO, NF90_CLOBBER)
        mode_flag = IOR(mode_flag, NF90_NETCDF4)
        status = nf90_create_par(rundata%filename, mode_flag, MPI_COMM_WORLD, info, file_id)
        if (status /= NF90_NOERR) then
            print *,'Could not open file ', rundata%filename
            return
        endif

        ! define the dimensions

        status = nf90_def_dim(file_id, 'X', rundata%globalnx, xdim_id)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' X'
        status = nf90_def_dim(file_id, 'Y', rundata%globalny, ydim_id)
        status = nf90_def_dim(file_id, 'velocity components', 2, vcomp_id)

        ! now that the dimensions are defined, we can define variables on them,...

        status = nf90_def_var(file_id, 'X coordinate', NF90_REAL, (/ xdim_id /), xcoord_id)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' xcoord'
        status = nf90_def_var(file_id, 'Y coordinate', NF90_REAL, (/ ydim_id /), ycoord_id)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' ycoord'

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

        ! To write out the values, we'll be using collective operations
        ! (NF90_INDEPENDANT is the other option, which can be better in some 
        ! situations)

        status = nf90_var_par_access(file_id, dens_id, NF90_COLLECTIVE)
        status = nf90_var_par_access(file_id, vel_id,  NF90_COLLECTIVE)

        allocate(x(rundata%globalnx), y(rundata%globalny))
        x = (/ (i-rundata%globalnx/2., i=1,rundata%globalnx) /)
        y = (/ (i-rundata%globalny/2., i=1,rundata%globalny) /)
        
        status = nf90_put_var(file_id, xcoord_id, x)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' X coord'
        status = nf90_put_var(file_id, ycoord_id, y)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Y coord'


        !
        ! Now we have to figure out the region within the global
        ! data that corresponds to our local data.
        !
        ! The n90_put_var() cals takes, to specify a region of an
        ! array, a vector of starting indicies, and a vector of counts
        ! drawing out the subregion.
        !
        !       |startx-|
        !       +-------|----|-------+   -+-
        !       |                    |    |
        !       |                    |  starty
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
        !  In this case the counts are (localnx,localny) and the offsets are 
        !  (startx,starty) = ((myx)/nxp*globalnx, (myy/nyp)*globalny)
        !

        densstarts(1) = (rundata % globalnx / rundata % npx) * rundata % myx + 1
        densstarts(2) = (rundata % globalny / rundata % npy) * rundata % myy + 1
        denscounts(1) = rundata % localnx 
        denscounts(2) = rundata % localny 

        print '(A,I3,A,I3,X,I3,X,I3,X,I3)', '[',rundata%rank,']: denstarts, denscounts = ', densstarts, denscounts
        status = nf90_put_var(file_id, dens_id, dens, start=densstarts, count=denscounts)
        if (status /= NF90_NOERR) print *, trim(nf90_strerror(status)), ' Dens'

        ! With the velocity, we have the complication of a third dimension, but
        ! it's a fairly simple one - we all have the same region (both velocity
        ! components)

        velstarts(1) = 1 
        velstarts(2) = (rundata % globalnx / rundata % npx) * rundata % myx + 1
        velstarts(3) = (rundata % globalny / rundata % npy) * rundata % myy  + 1
        velcounts(1) = 2
        velcounts(2) = rundata % localnx 
        velcounts(3) = rundata % localny 

        print '(A,I3,A,3(I3,X),3(I3,X))', '[',rundata%rank,']: velstarts, velcounts = ', velstarts, velcounts
        status = nf90_put_var(file_id, vel_id,  vel, start=velstarts, count=velcounts)
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
        i = 1

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

        sigma = gnx/4.

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
        double precision, intent(in), dimension(:,:) :: dens
        integer, intent(in) :: rank

        print *, 'Rank: [',rank,']', dens
    end subroutine printarray2d
  
    subroutine nearsquare(nprocs, npx, npy)
        implicit none
        integer, intent(IN)  :: nprocs
        integer, intent(OUT) :: npx, npy

        integer :: sq, n
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
