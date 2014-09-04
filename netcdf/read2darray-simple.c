#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include "netcdf.h"

#define MAXFILENAME 128
typedef struct rundata {
    size_t nx, ny, nveldims;
    char filename[MAXFILENAME]; 
} rundata_t;


double **array2d(int nx, int ny);
double ***array3d(int nd, int nx, int ny);

void readnetcdffile(rundata_t *rundata, double ***dens, double ****vel) {
    /* identifiers */
    int file_id, xdim_id, ydim_id, vcomp_id;
    int dens_id, vel_id;

    /* return status */
    int status;

    *dens = NULL;
    *vel  = NULL;

    /* Open a new file - read only */
    status = nc_open(rundata->filename, NC_NOWRITE, &file_id);

    /* netCDF routines return NC_NOERR on success */
    if (status != NC_NOERR) {
        fprintf(stderr,"Could not open file %s\n", rundata->filename);
        return;
    }

    /* Get the dimensions */
    status = nc_inq_dimid(file_id, "X", &xdim_id); 
    if (status != NC_NOERR) fprintf(stderr, "Could not get X\n");
    status = nc_inq_dimid(file_id, "Y", &ydim_id);
    if (status != NC_NOERR) fprintf(stderr, "Could not get Y\n");
    status = nc_inq_dimid(file_id, "velocity component", &vcomp_id);
    if (status != NC_NOERR) fprintf(stderr, "Could not get vcomp\n");

    status = nc_inq_dimlen(file_id, xdim_id, &(rundata->nx));
    if (status != NC_NOERR) fprintf(stderr, "Could not get X length\n");
    status = nc_inq_dimlen(file_id, ydim_id, &(rundata->ny));
    if (status != NC_NOERR) fprintf(stderr, "Could not get Y length\n");
    status = nc_inq_dimlen(file_id, vcomp_id, &(rundata->nveldims));
    if (status != NC_NOERR) fprintf(stderr, "Could not get vcomp length\n");

    /* now allocate the variables */
    *dens = array2d(rundata->nx, rundata->ny);
    *vel  = array3d(rundata->nveldims, rundata->nx, rundata->ny);

    nc_inq_varid(file_id, "Density", &dens_id);
    nc_inq_varid(file_id, "Velocity", &vel_id);
    
    /* Now read 'em in */
    nc_get_var_double(file_id, dens_id, &((*dens)[0][0]));
    nc_get_var_double(file_id, vel_id, &((*vel)[0][0][0]));

    nc_close(file_id);
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

    strcpy(rundata.filename,"data-simple.nc");
    get_options(argc,argv,&rundata);

    readnetcdffile(&rundata, &dens, &vel);
    printf("Read file %s: Density[%d][%d], Velocity[%d][%d][%d]\n", 
            rundata.filename, (int)rundata.nx, (int)rundata.ny,
            (int)rundata.nveldims, (int)rundata.nx, (int)rundata.ny);

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
                  puts("    --filename=S  (-f S): Set the output filename.");
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
