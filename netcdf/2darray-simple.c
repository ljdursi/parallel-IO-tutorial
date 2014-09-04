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
    int nx, ny;
    char filename[MAXFILENAME]; 
} rundata_t;


void writenetcdffile(rundata_t rundata, double **dens, double ***vel) {
    /* identifiers */
    int file_id, xdim_id, ydim_id, vcomp_id;
    int dens_id, vel_id;

    /* sizes */
    int densdims[2], veldims[3];

    /* return status */
    int status;

    /* Create a new file - clobber anything existing */
    status = nc_create(rundata.filename, NC_CLOBBER, &file_id);
    /* netCDF routines return NC_NOERR on success */
    if (status != NC_NOERR) {
        fprintf(stderr,"Could not open file %s\n", rundata.filename);
        return;
    }

    /* define the dimensions */
    nc_def_dim(file_id, "X", rundata.nx, &xdim_id);
    nc_def_dim(file_id, "Y", rundata.ny, &ydim_id);
    nc_def_dim(file_id, "velocity component", 2, &vcomp_id);

    /* now that we've defined the dimensions, we can define variables on them */
    densdims[0] = xdim_id;  densdims[1] = ydim_id;
    veldims[0]  = vcomp_id; veldims[1] = xdim_id; veldims[2] = ydim_id;

    nc_def_var(file_id, "Density",  NC_DOUBLE, 2, densdims, &dens_id);
    nc_def_var(file_id, "Velocity", NC_DOUBLE, 3, veldims,  &vel_id);

    /* we are now done defining variables and their attributes */
    nc_enddef(file_id);

    /* Write out the data to the variables we've defined */
    nc_put_var_double(file_id, dens_id, &(dens[0][0]));
    nc_put_var_double(file_id, vel_id,  &(vel[0][0][0]));

    
    nc_close(file_id);
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
    strcpy(rundata.filename,"data-simple.nc");

    get_options(argc,argv,&rundata);

    dens = array2d(rundata.nx, rundata.ny);
    vel  = array3d(2, rundata.nx, rundata.ny);
 
    fillarray2d(dens, rundata.nx, rundata.ny);
    fillarray3d(vel, rundata.nx, rundata.ny);

    if (rundata.nx*rundata.ny < 200) 
        printarray2d(dens, rundata.nx, rundata.ny);
    writenetcdffile(rundata, dens, vel);

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
                      if (!(tst=fopen(rundata->filename,"w"))) {
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
                  puts("    --filename=S   (-f S): Set the output filename.");
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
