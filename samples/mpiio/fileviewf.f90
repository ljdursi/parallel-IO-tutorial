program fileview
    use mpi
    implicit none

    integer(mpi_offset_kind) :: offset
    integer, dimension(mpi_status_size) :: wstatus
    integer :: fileno
    integer, parameter :: datasize=15
    integer :: ierr, rank, comsize
    integer, dimension(2) :: gridsize
    integer ::  myrow, mycol, locnrows, locncols
    integer ::  startr, endr, startc, endc
    character, allocatable, dimension (:,:) :: mydata
    integer :: viewtype

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    call MPI_Dims_create(comsize, 2, gridsize, ierr)
    mycol = mod(rank, gridsize(2))
    myrow = rank/gridsize(2)
  
    ! start, end, start at 0 even for Fortran
    locncols = (datasize+1) / gridsize(2)
    startc = mycol * locncols
    endc   = startc + locncols - 1
    if (mycol == gridsize(2) - 1) then
        endc = datasize
        locncols = endc - startc + 1
    endif
  
    ! start, end, start at 0 even for Fortran
    locnrows = datasize / gridsize(1)
    startr = myrow * locnrows
    endr   = startr + locnrows - 1
    if (myrow == gridsize(1) - 1) then
        endr = datasize - 1
        locnrows = endr - startr + 1
    endif
  
    print '(A,I3,A,2(I3),A,2(I3),2(I3))', 'Rank: ', rank, ' size: ', locnrows,  locncols, ' starts: ', startr, endr, startc, endc
    print '(3(A,I3))', 'Rank: ', rank, ' myrow: ', myrow,  ' mycol: ', mycol

    allocate(mydata(locncols,locnrows))
    mydata = achar(ichar('0') + rank)
    if (mycol == gridsize(2) - 1) mydata(locncols,:) = achar(10)
    
    call MPI_Type_create_subarray(2, [datasize, datasize+1], [locnrows, locncols], &
                                 [startr, startc], MPI_ORDER_C, &
                                 MPI_CHARACTER, viewtype, ierr)
    call MPI_Type_commit(viewtype, ierr)

    offset = 0

    call MPI_File_open(MPI_COMM_WORLD, "viewtype.txt",&
                       ior(MPI_MODE_CREATE,MPI_MODE_WRONLY),   &
                       MPI_INFO_NULL, fileno, ierr)
    call MPI_File_set_view(fileno, offset,  MPI_CHARACTER, viewtype, &
                           "native", MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fileno, mydata, locnrows*locncols, &
                        MPI_CHARACTER, wstatus, ierr)

    deallocate(mydata)
    call MPI_Type_free(viewtype, ierr)

    call MPI_File_close(fileno, ierr)
   
    call MPI_Finalize(ierr)

end program fileview
