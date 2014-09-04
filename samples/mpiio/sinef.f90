program sinef
    use mpi
    implicit none

    integer(mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: wstatus
    integer :: fileno
    integer, parameter :: npts=200
    integer :: ierr, rank, comsize
    integer ::  start, locnpts
    real, allocatable, dimension (:) :: mydata
    integer :: viewtype
    integer :: i

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! start, end, start at 0 even for Fortran
    locnpts = npts/comsize
    start   = locnpts*rank
    if (rank == comsize-1) then
        locnpts = npts - start
    endif
  
    allocate(mydata(locnpts))
    do i=1,locnpts 
        mydata(i) = sin((start+i)*1.0*8*atan(1.)/npts)
    enddo 
    
    call MPI_File_open(MPI_COMM_WORLD, "sine.dat",              &
                       ior(MPI_MODE_CREATE,MPI_MODE_WRONLY),    &
                       MPI_INFO_NULL, fileno, ierr)
    
    !! other mpi-io stuff

    call MPI_File_close(fileno, ierr)
   
    call MPI_Finalize(ierr)

end program sinef
