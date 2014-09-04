#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "hdf5.h"


#define MAXFILENAME 128
typedef struct rundata {
    int nx, ny, nveldims;
    char filename[MAXFILENAME]; 
} rundata_t;

double **array2d(int nx, int ny);
double ***array3d(int nd, int nx, int ny);

void readhdf5file(rundata_t *rundata, double ***dens, double ****vel) {
    /* identifiers */
    hid_t file_id, dens_dataset_id, vel_dataset_id;
    hid_t dens_dataspace_id, vel_dataspace_id;

    /* sizes */
    hsize_t densdims[2], veldims[3];

    /* status */
    herr_t status;

    printf("Attempting to open file %s\n", rundata->filename);
    /* Open a file read only, with default file access properties */
    file_id = H5Fopen(rundata->filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* HDF5 routines generally return a negative number on failure.  
     * Should check return values! */
    if (file_id < 0) {
        fprintf(stderr,"Could not open file %s\n", rundata->filename);
        return;
    }

    /* Open the datasets from the file. */
    dens_dataset_id = H5Dopen(file_id, "dens", H5P_DEFAULT);
    vel_dataset_id  = H5Dopen(file_id, "vel",  H5P_DEFAULT);
    if (dens_dataset_id < 0 || vel_dataset_id < 0) {
        fprintf(stderr,"Could not open datasets!\n");
        return;
    }

    /* get dimensions */
    dens_dataspace_id = H5Dget_space(dens_dataset_id);
    if (dens_dataspace_id < 0) {
        fprintf(stderr,"Could not get dens dataspace id!\n");
        return;
    }
    status = H5Sget_simple_extent_dims(dens_dataspace_id, densdims, NULL);
    if (status != 2) {
        fprintf(stderr,"Could not get dens extents!\n");
        return;
    }

    vel_dataspace_id = H5Dget_space(vel_dataset_id);
    if (vel_dataspace_id < 0) {
        fprintf(stderr,"Could not get vel dataspace id!\n");
        return;
    }
    status = H5Sget_simple_extent_dims(vel_dataspace_id, veldims, NULL);
    if (status != 3) {
        fprintf(stderr,"Could not get vel extents!\n");
        return;
    }

    rundata -> nx = densdims[0];
    rundata -> ny = densdims[1];
    rundata -> nveldims = veldims[0];

    /* Make data arrays */

    *dens = array2d(densdims[0], densdims[1]);
    *vel  = array3d(veldims[0],  veldims[1],  veldims[2]);

    /* Read the data.  */
    
    H5Dread(dens_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &((*dens)[0][0]));
    H5Dread(vel_dataset_id,  H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &((*vel)[0][0][0]));

    /* End access to groups & data sets and release resources used by them */
    status = H5Sclose(dens_dataspace_id);
    status = H5Dclose(dens_dataset_id);
    status = H5Sclose(vel_dataspace_id);
    status = H5Dclose(vel_dataset_id);
    
    /* Close the file */
    status = H5Fclose(file_id);
    return;
}



int get_options(int argc, char **argv, rundata_t *rundata);
void freearray2d(double **d);
void freearray3d(double ***d);
void fillarray2d(double **d, int nx, int ny);
void fillarray3d(double ***d, int nx, int ny);
void printarray2d(double **d, int nx, int ny);


int main(int argc, char **argv) {
    double **dens;
    double ***vel;
    rundata_t rundata;

    strcpy(rundata.filename,"data-simple.h5");
    get_options(argc,argv,&rundata);

    readhdf5file(&rundata, &dens, &vel);
    printf("Read file %s: Density[%d][%d], Velocity[%d][%d][%d]\n", 
            rundata.filename, rundata.nx, rundata.ny,
            rundata.nveldims, rundata.nx, rundata.ny);
   
    if (rundata.nx*rundata.ny < 200) 
        printarray2d(dens, rundata.nx, rundata.ny);

    freearray2d(dens);
    freearray3d(vel);

    return 0;
}  


int get_options(int argc, char **argv, rundata_t *rundata) {
                                        
    struct option long_options[] = {
        {"nx",       required_argument, 0, 'x'},
        {"ny",       required_argument, 0, 'y'},
        {"filename", required_argument, 0, 'f'},
        {"help",     no_argument, 0, 'h'},
        {0, 0, 0, 0}};

    char c;
    int option_index;
    FILE *tst;
    char *defaultfname=rundata->filename;

    while (1) { 
        c = getopt_long(argc, argv, "x:y:f:h", long_options, 
                        &option_index);
        if (c == -1) break;

        switch (c) { 
            case 0: if (long_options[option_index].flag != 0)
                    break;

            case 'f': strncpy(rundata->filename, optarg, MAXFILENAME-1);
                      if (!(tst=fopen(rundata->filename,"r"))) {
                          fprintf(stderr,
                                  "Cannot use filename %s;\n",
                                  rundata->filename);
                          fprintf(stderr, "  Using %s\n",defaultfname);
                          strcpy(rundata->filename, defaultfname);
                      } else 
                          fclose(tst);
                      break; 

            case 'h':
                  puts("Options: ");
                  puts("    --fileaname=S  (-f S): Set the filename.");
                  puts("");
                  return +1;

            default: printf("Invalid option %s\n", optarg);
                 break;
        }
    }
    return 0;
}

double **array2d(int nx, int ny) {
    int i;
    double *data = (double *)malloc(nx*ny*sizeof(double));
    double **p = (double **)malloc(nx*sizeof(double *));

    if (data == NULL) return NULL;
    if (p == NULL) {
        free(data);
        return NULL;
    }

    for (i=0; i<nx; i++) {
        p[i] = &(data[ny*i]);
    }
    
    return p;
}

void freearray2d(double **p) {
    free(p[0]);
    free(p);
    return;
}

double ***array3d(int nd, int nx, int ny) {
    int i;
    double *data = (double *)malloc(nd*nx*ny*sizeof(double));
    double **datap = (double **)malloc(nd*nx*sizeof(double *));
    double ***p = (double ***)malloc(nd*sizeof(double **));

    if (data == NULL) return NULL;
    if (datap == NULL) {
        free(data);
        return NULL;
    }
    if (p == NULL) {
        free(data);
        free(datap);
        return NULL;
    }

    for (i=0; i<nd*nx; i++) {
        datap[i] = &(data[ny*i]);
    }

    for (i=0; i<nd; i++) {
        p[i] = &(datap[nx*i]);
    }
    
    return p;
}

void freearray3d(double ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
    return;
}

void printarray2d(double **d, int nx, int ny) {
    int i,j;

    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            printf("%10.3g\t",d[i][j]);
        }
        puts("");
    } 

    return;
}
