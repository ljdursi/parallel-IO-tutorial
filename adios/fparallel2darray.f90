program parallelarray
    use mpi
    implicit none

    type :: rundata_t
        integer :: globalnx, globalny
        integer :: localnx, localny
        integer :: npx, npy
        integer :: myx, myy
        integer :: startx, starty
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
    call writeadiosfile(rundata, dens, vel)

    call MPI_Finalize(ierr)

    deallocate(dens)
    deallocate(vel)


contains
    subroutine writeadiosfile(rundata, dens, vel)

        use iso_fortran_env
        implicit none
        type(rundata_t), intent(in) :: rundata
        double precision, intent(in), dimension(:,:) :: dens
        double precision, intent(in), dimension(:,:,:) :: vel

        integer :: adios_err
        integer(kind=int64) :: adios_groupsize
        integer(kind=int64) :: adios_totalsize
        integer(kind=int64) :: adios_handle
        integer :: comm, comsize, ierr
 
        comm = MPI_COMM_WORLD
        call MPI_Comm_size(comm, comsize, ierr)

        call adios_init ("fadios_global.xml", adios_err)
        call adios_open (adios_handle, "ArrayData", rundata%filename, "w", comm, adios_err)
        include "gwrite_ArrayData.fh"
        call adios_close (adios_handle, adios_err)

    end subroutine writeadiosfile

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
 
        rundata % startx = ((rundata % globalnx)/(rundata % npx))*(rundata % myx)
        rundata % starty = ((rundata % globalny)/(rundata % npy))*(rundata % myy)

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

        startx = rundata % startx
        starty = rundata % starty

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

        startx = rundata % startx
        starty = rundata % starty

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
