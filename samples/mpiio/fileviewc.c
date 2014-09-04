#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

char **chararray2d_alloc(int n, int m) {
    char *data = malloc(n*m*sizeof(char));
    char **ptrs = malloc(n*sizeof(char*));
    int i;
    for (i=0; i<n; i++) 
        ptrs[i] = &(data[i*m]);

    return ptrs;
}

void chararray2d_free(int n, int m, char **a) {
    free(&(a[0][0]));
    free(a);
}

int main(int argc, char **argv) {
    int ierr, rank, size;
    MPI_Offset offset;
    MPI_File   file;
    MPI_Status status;
    const int datasize=15;
    int gridsize[2];
    int globalsize[2] = {datasize, datasize+1};
    int subsize[2], start[2];
    int  myrow, mycol, locnrows, locncols;
    int startr, endr, startc, endc;
    char **mydata;
    MPI_Datatype viewtype;
    int i, j;

    ierr = MPI_Init(&argc, &argv);
    ierr|= MPI_Comm_size(MPI_COMM_WORLD, &size);
    ierr|= MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    gridsize[0]=0; gridsize[1] = 0;
    ierr = MPI_Dims_create(size, 2, gridsize);
    
    printf("Rank = %d, globalsize = (%d, %d), gridsize = (%d, %d)\n", rank, globalsize[0], globalsize[1], gridsize[0], gridsize[1]);

    mycol = rank % gridsize[1];
    myrow = rank / gridsize[1];
  
    locncols = globalsize[1] / gridsize[1];
    startc = mycol * locncols;
    endc   = startc + locncols - 1;
    if (mycol == gridsize[1] - 1) {
        endc = datasize;
        locncols = endc - startc + 1;
    }
  
    locnrows = globalsize[0] / gridsize[0];
    startr = myrow * locnrows;
    endr   = startr + locnrows - 1;
    if (myrow == gridsize[0] - 1) {
        endr = datasize - 1;
        locnrows = endr - startr + 1;
    }
  
    printf("Rank = %d, size = (%d, %d), starts = (%d, %d, %d, %d)\n", rank, locnrows, locncols, startr, endr, startc, endc);

    mydata = chararray2d_alloc(locnrows, locncols);
    for (i=0; i<locnrows; i++)
        for (j=0; j<locncols; j++) 
            mydata[i][j] = (char)('0' + rank);
    
    if (mycol == gridsize[1] - 1) 
        for (i=0; i<locnrows; i++) 
            mydata[i][locncols-1] = '\n';

    subsize[0] = locnrows; subsize[1] = locncols;
    start[0]   = startr;   start[1]   = startc;

    MPI_Type_create_subarray(2, globalsize, subsize, start, MPI_ORDER_C,
                                 MPI_CHAR, &viewtype);
    MPI_Type_commit(&viewtype);

    offset = 0;

    MPI_File_open(MPI_COMM_WORLD, "viewtype.txt", 
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    MPI_File_set_view(file, offset,  MPI_CHAR, viewtype, 
                           "native", MPI_INFO_NULL);

    MPI_File_write_all(file, &(mydata[0][0]), locnrows*locncols, MPI_CHAR, &status);
    MPI_File_close(&file);
   
    MPI_Type_free(&viewtype);
    chararray2d_free(locnrows, locncols, mydata);

    MPI_Finalize();
    return 0;
}
