program noncontig
    use mpi
    implicit none

    integer(mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: wstatus
    integer :: fileno
    integer, parameter :: msgsize=6, strsize=12
    character(strsize) :: message
    integer :: ierr, rank, comsize
    integer :: everyother

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    if (mod(rank,2) == 0) then
        message = "H@e#l*l^o* A"
    else 
        message = "WFoQr#l>d@!_"
    endif

    offset = rank*msgsize

    call MPI_Type_vector(msgsize, 1, 2, MPI_CHARACTER, everyother, ierr)
    call MPI_Type_commit(everyother, ierr)

    call MPI_File_open(MPI_COMM_WORLD, "helloworld-nc.txt",     &
                       ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
                       MPI_INFO_NULL, fileno, ierr)

    call MPI_File_seek (fileno, offset, MPI_SEEK_SET, ierr)
    call MPI_File_write(fileno, message, 1, everyother, &
                        wstatus, ierr)
    call MPI_File_close(fileno, ierr)
   
    call MPI_Finalize(ierr)

end program noncontig
