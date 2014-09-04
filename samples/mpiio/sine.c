#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv) {
    int ierr, rank, size;
    MPI_Offset offset;
    MPI_File   file;
    MPI_Status status;
    int npts=200;
    int start;
    int locnpts;
    float *data;
    MPI_Datatype viewtype;

    ierr = MPI_Init(&argc, &argv);
    ierr|= MPI_Comm_size(MPI_COMM_WORLD, &size);
    ierr|= MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    locnpts = npts/size;
    start = locnpts * rank;
    if (rank == size-1)
        locnpts = (npts-start);

    data = malloc(locnpts * sizeof(float));

    for (int i=0; i<locnpts; i++)
        data[i] = sin((start+i)*1.0*8*atan(1.)/npts);

    MPI_File_open(MPI_COMM_WORLD, "sine.dat", 
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    /*  Other MPI-IO stuff */

    MPI_File_write_all(file, data, locnpts, MPI_FLOAT, &status);

    MPI_File_close(&file);
   
    free(data);
    MPI_Finalize();
    return 0;
}
