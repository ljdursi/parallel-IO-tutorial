program MPIIO_helloworld
    use mpi
    implicit none

    integer(mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: wstatus
    integer :: fileno
    integer, parameter :: msgsize=6
    character(msgsize) :: message
    integer :: ierr, rank, comsize

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    if (mod(rank,2) == 0) then
        message = "Hello "
    else 
        message = "World!"
    endif

    offset = rank*msgsize

    call MPI_File_open(MPI_COMM_WORLD, "helloworld-at-all.txt",&
                       ior(MPI_MODE_CREATE,MPI_MODE_WRONLY),   &
                       MPI_INFO_NULL, fileno, ierr)

    call MPI_File_write_at_all(fileno, offset, message, msgsize, &
                        MPI_CHARACTER, wstatus, ierr)
    call MPI_File_close(fileno, ierr)
   
    call MPI_Finalize(ierr)

end program MPIIO_helloworld
