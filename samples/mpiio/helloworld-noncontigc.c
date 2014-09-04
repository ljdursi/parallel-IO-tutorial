#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv) {
    int ierr, rank, size;
    MPI_Offset offset;
    MPI_File   file;
    MPI_Status status;
    MPI_Datatype everyother;
    const int msgsize=6;
    const int strsize=12;
    char message[strsize+1];

    ierr = MPI_Init(&argc, &argv);
    ierr|= MPI_Comm_size(MPI_COMM_WORLD, &size);
    ierr|= MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if ((rank % 2) == 0)
        strcpy (message, "H@e#l&l^o* A"); 
    else    
        strcpy (message, "WFoQr#l>d@!_");

    printf("Rank %d has message <%s>\n", rank, message);

    offset = (msgsize*rank);

    MPI_Type_vector(msgsize, 1, 2, MPI_CHAR, &everyother);
    MPI_Type_commit(&everyother);

    MPI_File_open(MPI_COMM_WORLD, "helloworld-nc.txt", 
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);
    MPI_File_seek(file, offset, MPI_SEEK_SET);
    MPI_File_write(file, message, 1, everyother, &status);
    MPI_File_close(&file);

    MPI_Type_free(&everyother);
   
    MPI_Finalize();
    return 0;
}

