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
    int nx, ny;
    char filename[MAXFILENAME]; 
} rundata_t;


void writehdf5file(rundata_t rundata, double **dens, double ***vel) {
    /* identifiers */
    hid_t file_id, arr_group_id, dens_dataset_id, vel_dataset_id;
    hid_t dens_dataspace_id, vel_dataspace_id;

    /* sizes */
    hsize_t densdims[2], veldims[3];

    /* status */
    herr_t status;

    /* Create a new file - truncate anything existing, use default properties */
    file_id = H5Fcreate(rundata.filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* HDF5 routines generally return a negative number on failure.  
     * Should check return values! */
    if (file_id < 0) {
        fprintf(stderr,"Could not open file %s\n", rundata.filename);
        return;
    }

    /* Create a new group within the new file */
    arr_group_id = H5Gcreate(file_id,"/ArrayData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Give this group an attribute listing the time of calculation */
    {
        hid_t attr_id,attr_sp_id;
        struct tm *t;
        time_t now;
        int yyyymm;
        now = time(NULL);
        t = localtime(&now);
        yyyymm = (1900+t->tm_year)*100+t->tm_mon;

        attr_sp_id = H5Screate(H5S_SCALAR);
        attr_id = H5Acreate(arr_group_id, "Calculated on (YYYYMM)", H5T_STD_U32LE, attr_sp_id, H5P_DEFAULT, H5P_DEFAULT);
        printf("yymm = %d\n",yyyymm);
        H5Awrite(attr_id, H5T_NATIVE_INT, &yyyymm);
        H5Aclose(attr_id);
        H5Sclose(attr_sp_id);
    }


    /* Create the data space for the two datasets. */
    densdims[0] = rundata.nx; densdims[1] = rundata.ny;
    veldims[0] = 2; veldims[1] = rundata.nx; veldims[2] = rundata.ny;
    
    dens_dataspace_id = H5Screate_simple(2, densdims, NULL);
    vel_dataspace_id  = H5Screate_simple(3, veldims,  NULL);

    /* Create the datasets within the file. 
     * H5T_IEEE_F64LE is a standard (IEEE) double precision (64 bit) floating (F) data type
     * and will work on any machine.  H5T_NATIVE_DOUBLE would work too, but woudl give
     * different results on GPC and TCS */

    dens_dataset_id = H5Dcreate(file_id, "/ArrayData/dens", H5T_IEEE_F64LE, 
                                dens_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    vel_dataset_id  = H5Dcreate(file_id, "/ArrayData/vel",  H5T_IEEE_F64LE, 
                                vel_dataspace_id,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Write the data.  We're writing it from memory, where it is saved 
     * in NATIVE_DOUBLE format */
    status = H5Dwrite(dens_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(dens[0][0]));
    status = H5Dwrite(vel_dataset_id,  H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vel[0][0][0]));

    /* We'll create another group for related info and put some things in there */
    {   
        hid_t other_group_id;
        hid_t timestep_id, timestep_space;
        hid_t comptime_id, comptime_space;
        hid_t author_id, author_space, author_type;
        char *authorname="Jonathan Dursi";
        int timestep=13;
        float comptime=81.773;

        /* create group */
        other_group_id = H5Gcreate(file_id,"/OtherStuff", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        /* scalar space, data for integer timestep */
        timestep_space = H5Screate(H5S_SCALAR);
        timestep_id = H5Dcreate(other_group_id, "Timestep", H5T_STD_U32LE,
                                timestep_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(timestep_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &timestep);
        H5Dclose(timestep_id);
        H5Sclose(timestep_space);

        /* scalar space, data for floating compute time */
        comptime_space = H5Screate(H5S_SCALAR);
        comptime_id = H5Dcreate(other_group_id, "Compute Time", H5T_IEEE_F32LE,
                                comptime_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(comptime_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &comptime);
        H5Dclose(comptime_id);
        H5Sclose(comptime_space);

        /* scalar space, data for author name */
        author_space = H5Screate(H5S_SCALAR);
        author_type  = H5Tcopy(H5T_C_S1);   /* copy the character type.. */
        status = H5Tset_size (author_type, strlen(authorname));  /* and make it longer */
        author_id = H5Dcreate(other_group_id, "Simulator Name", author_type, author_space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(author_id, author_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, authorname);
        H5Dclose(author_id);
        H5Sclose(author_space);
        H5Tclose(author_type);

        H5Gclose(other_group_id);
    }

    /* End access to groups & data sets and release resources used by them */
    status = H5Sclose(dens_dataspace_id);
    status = H5Dclose(dens_dataset_id);
    status = H5Sclose(vel_dataspace_id);
    status = H5Dclose(vel_dataset_id);
    status = H5Gclose(arr_group_id);
    
    /* Close the file */
    status = H5Fclose(file_id);
    return;
}



int get_options(int argc, char **argv, rundata_t *rundata);
double **array2d(int nx, int ny);
double ***array3d(int nd, int nx, int ny);
void freearray2d(double **d);
void freearray3d(double ***d);
void fillarray2d(double **d, int nx, int ny);
void fillarray3d(double ***d, int nx, int ny);
void printarray2d(double **d, int nx, int ny);


int main(int argc, char **argv) {
    double **dens;
    double ***vel;
    rundata_t rundata;

    rundata.nx = 100;
    rundata.ny = 100;
    strcpy(rundata.filename,"data.h5");

    get_options(argc,argv,&rundata);

    dens = array2d(rundata.nx, rundata.ny);
    vel  = array3d(2, rundata.nx, rundata.ny);
 
    fillarray2d(dens, rundata.nx, rundata.ny);
    fillarray3d(vel, rundata.nx, rundata.ny);

    if (rundata.nx*rundata.ny < 200) 
        printarray2d(dens, rundata.nx, rundata.ny);
    writehdf5file(rundata, dens, vel);

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
    int tempint;
    int defaultnpts = 100;
    FILE *tst;
    char *defaultfname=rundata->filename;

    while (1) { 
        c = getopt_long(argc, argv, "x:y:f:h", long_options, 
                        &option_index);
        if (c == -1) break;

        switch (c) { 
            case 0: if (long_options[option_index].flag != 0)
                    break;

            case 'x': tempint = atoi(optarg);
                      if (tempint < 1 || tempint > 500) {
                          fprintf(stderr,
                                  "%s: Cannot use number of points %s;\n",
                                  argv[0], optarg);
                          fprintf(stderr,"  Using %d\n", defaultnpts);
                          rundata->nx = defaultnpts;
                      } else {
                          rundata->nx = tempint;
                      }
                      break;

            case 'y': tempint = atoi(optarg);
                      if (tempint < 1 || tempint > 500) {
                          fprintf(stderr,
                                  "%s: Cannot use number of points %s;\n",
                                  argv[0], optarg);
                          fprintf(stderr,"  Using %d\n", defaultnpts);
                          rundata->ny = defaultnpts;
                      } else {
                          rundata->ny = tempint;
                      }
                      break;

            case 'f': strncpy(rundata->filename, optarg, MAXFILENAME-1);
                      if (!(tst = fopen(rundata->filename,"w"))) {
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
                  puts("    --nx=N         (-x N): Set the number of grid cells in x direction.");
                  puts("    --ny=N         (-y N): Set the number of grid cells in y direction.");
                  puts("    --fileaname=S  (-f S): Set the output filename.");
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

void fillarray2d(double **d, int nx, int ny) {
    int i,j;
    double sigma=nx/4.;
    double r2,r2x;

    for (i=0;i<nx;i++) {
        r2x = (i-nx/2.)*(i-nx/2.);
        for (j=0;j<ny;j++) {
            r2 = r2x + (j-ny/2.)*(j-ny/2.);
            d[i][j] = 1. + 4.*exp(-r2/(2.*sigma*sigma));
        }
    } 

    return;
}

void fillarray3d(double ***d, int nx, int ny) {
    int i,j;
    double sigma=nx/4.;
    double r2,r2x;

    for (i=0;i<nx;i++) {
        r2x = (i-nx/2.)*(i-nx/2.);
        for (j=0;j<ny;j++) {
            r2 = r2x + (j-ny/2.)*(j-ny/2.);
            if (r2 == 0) {
                d[0][i][j] = 0.;
                d[1][i][j] = 0.;
            } else {
                d[0][i][j] = exp(-r2/(2.*sigma*sigma))*(j-ny/2.);
                d[1][i][j] = -exp(-r2/(2.*sigma*sigma))*(i-nx/2.);
            } 
        }
    } 

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

