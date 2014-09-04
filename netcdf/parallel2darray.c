#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "netcdf.h"

#define MAXFILENAME 128
typedef struct rundata {
    int globalnx, globalny;
    int localnx, localny;
    int npx, npy, nprocs, rank;
    int myx, myy;
    char filename[MAXFILENAME]; 
} rundata_t;

void writenetcdffile(rundata_t rundata, double **dens, double ***vel) {
    /* identifiers */
    int file_id, xdim_id, xcoord_id, ydim_id, ycoord_id, vcomp_id;
    int dens_id, vel_id;

    /* sizes */
    int densdims[2], veldims[3];

    /* array containing x, y coordinates */
    float *x, *y;
    int i;
    const char *coordunit="cm";

    /* name of units for dens, vel */
    const char *densunit="g/cm^3";
    const char *velunit ="cm/s";

    /* offsets for sub-regions of arrays */
    size_t starts[3];
    size_t counts[3];

    /* return status */
    int status;

    /* MPI-IO hints for performance */
    MPI_Info info;

    /* set up x, y coordinates */
    x = (float *)malloc(rundata.globalnx * sizeof(float));
    y = (float *)malloc(rundata.globalny * sizeof(float));
    for (i=0; i<rundata.globalnx; i++) 
        x[i] = (1.*i-rundata.globalnx/2.);
    for (i=0; i<rundata.globalny; i++) 
        y[i] = (1.*i-rundata.globalny/2.);
    
    
    /* set the MPI-IO hints for better performance on GPFS */
    MPI_Info_create(&info);
    MPI_Info_set(info,"IBM_largeblock_io","true");

    /* Create a new file - clobber anything existing */
    status = nc_create_par(rundata.filename, NC_MPIIO | NC_CLOBBER | NC_NETCDF4, 
                       MPI_COMM_WORLD, info, &file_id);
    /* netCDF routines return NC_NOERR on success */
    if (status != NC_NOERR) {
        fprintf(stderr,"Could not open file %s\n", rundata.filename);
        return;
    }

    /* define the dimensions */
    nc_def_dim(file_id, "X", rundata.globalnx, &xdim_id);
    nc_def_dim(file_id, "Y", rundata.globalny, &ydim_id);
    nc_def_dim(file_id, "velocity component", 2, &vcomp_id);

    /* define the coordinate variables,... */
    nc_def_var(file_id, "X coordinate", NC_FLOAT, 1, &xdim_id, &xcoord_id);
    nc_def_var(file_id, "Y coordinate", NC_FLOAT, 1, &ydim_id, &ycoord_id);

    /* ...and assign units to them as an attribute */
    nc_put_att_text(file_id, xcoord_id, "units", strlen(coordunit), coordunit);
    nc_put_att_text(file_id, ycoord_id, "units", strlen(coordunit), coordunit);

    /* now that we've defined the dimensions, we can define variables on them */
    densdims[0] = xdim_id;  densdims[1] = ydim_id;
    veldims[0]  = vcomp_id; veldims[1] = xdim_id; veldims[2] = ydim_id;

    nc_def_var(file_id, "Density",  NC_DOUBLE, 2, densdims, &dens_id);
    nc_def_var(file_id, "Velocity", NC_DOUBLE, 3, veldims,  &vel_id);

    /* assign units to the variables */
    nc_put_att_text(file_id, dens_id, "units", strlen(densunit), densunit);
    nc_put_att_text(file_id, vel_id,  "units", strlen(velunit),  velunit);

    /* we are now done defining variables and their attributes */
    nc_enddef(file_id);

    /* Write out the data to the variables we've defined */
    nc_put_var_float(file_id, xcoord_id, x);
    nc_put_var_float(file_id, ycoord_id, y);
    /*
     * Now we have to figure out the region within the global
     * data that corresponds to our local data.
     *
     * The nc_put_vara_type() takes, to specify a region of an
     * array, a vector of starting indicies, and a vector of counts
     * drawing out the subregion.
     *
     *       |startx-|
     *       +-------|----|-------+   -+-
     *       |                    |    |
     *       |                    |  starty
     *       |                    |    |
     *       -       +----+       -   -+-
     *       |       |    |       |    |
     *       |       |    |       |  localny
     *       |       |    |       |    |
     *       -       +----+       -   -+-
     *       |                    |
     *       |                    |
     *       +-------|----|-------+
     *               localnx
     *
     *  In this case the blocksizes are (localnx,localny) and the offsets are 
     *  (offx,offy) = ((myx)/nxp*globalnx, (myy/nyp)*globalny)
     */

    /* The big data will be written to collectively; 
     * the alternative is NC_INDEPENDENT  */
    nc_var_par_access(file_id, dens_id, NC_COLLECTIVE);
    nc_var_par_access(file_id, vel_id,  NC_COLLECTIVE);

    /* densities */
    starts[0] = (rundata.globalnx/rundata.npx)*rundata.myx;
    starts[1] = (rundata.globalny/rundata.npy)*rundata.myy;
    counts[0]  = rundata.localnx;
    counts[1]  = rundata.localny;

    nc_put_vara_double(file_id, dens_id, starts, counts, &(dens[0][0]));

    /* velocities */
    starts[0] = 0;
    starts[1] = (rundata.globalnx/rundata.npx)*rundata.myx;
    starts[2] = (rundata.globalny/rundata.npy)*rundata.myy;
    counts[0] = 2;
    counts[1] = rundata.localnx;
    counts[2] = rundata.localny;

    nc_put_vara_double(file_id, vel_id, starts, counts, &(vel[0][0][0]));

    nc_close(file_id);
    return;
}


int get_options(int argc, char **argv, rundata_t *rundata);
double **array2d(int nx, int ny);
double ***array3d(int nd, int nx, int ny);
void freearray2d(double **d);
void freearray3d(double ***d);
void fillarray2d(double **d, rundata_t *r);
void fillarray3d(double ***v, rundata_t *r);
void printarray2d(double **d, int nx, int ny);
void nearsquare(int nprocs, int *npx, int *npy);


int main(int argc, char **argv) {
    double **locdens;
    double ***locvel;
    rundata_t rundata;
    int ierr;
    int rank, size;


    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    /*  
     * set default values for parameters, then check what the user says
     * with the options 
     */

    /* choose sensible defaults */
    nearsquare(size,&(rundata.npx),&(rundata.npy));
    rundata.globalnx = 100;
    rundata.globalny = 100;
    rundata.nprocs = size;
    rundata.rank = rank;
    rundata.localnx = (rundata.globalnx / rundata.npx);
    rundata.localny = (rundata.globalny / rundata.npy);
    strcpy(rundata.filename,"paralleldata.nc");

    /* get options */
    get_options(argc,argv,&rundata);

    printf("[%d]: Assigned to region (%d,%d) \n", rank, rundata.myx, rundata.myy);
    
    /*
     * allocate our local arrays
     */
    locdens = array2d(rundata.localnx, rundata.localny);
    locvel  = array3d(2, rundata.localnx, rundata.localny);
    printf("[%d]: Allocated arrays\n", rank);
 
    fillarray2d(locdens, &rundata);
    fillarray3d(locvel, &rundata);
    printf("[%d]: Filled arrays\n", rank);

    if (rundata.localnx*rundata.localny < 200) 
        printarray2d(locdens, rundata.localnx, rundata.localny);
    writenetcdffile(rundata, locdens, locvel);
    printf("[%d]: Wrote file\n", rank);

    freearray2d(locdens);
    freearray3d(locvel);

    ierr = MPI_Finalize();
    return 0;
}  


int get_options(int argc, char **argv, rundata_t *rundata) {
                                        
    struct option long_options[] = {
        {"nx",       required_argument, 0, 'x'},
        {"ny",       required_argument, 0, 'y'},
        {"npx",      required_argument, 0, 'X'},
        {"npy",      required_argument, 0, 'Y'},
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
                          rundata->globalnx = defaultnpts;
                      } else {
                          rundata->globalnx = tempint;
                      }
                      break;

            case 'y': tempint = atoi(optarg);
                      if (tempint < 1 || tempint > 500) {
                          fprintf(stderr,
                                  "%s: Cannot use number of points %s;\n",
                                  argv[0], optarg);
                          fprintf(stderr,"  Using %d\n", defaultnpts);
                          rundata->globalny = defaultnpts;
                      } else {
                          rundata->globalny = tempint;
                      }
                      break;


            case 'X': tempint = atoi(optarg);
                      if (tempint < 1 || tempint > rundata->nprocs) {
                          fprintf(stderr,
                                  "%s: Cannot use number of processors in x direction %s;\n",
                                  argv[0], optarg);
                          fprintf(stderr,"  Using %d\n", rundata->npx);
                      } else if (rundata->nprocs % tempint != 0) {
                          fprintf(stderr,
                                  "%s: Number of processors in x direction %s does not divide %d;\n",
                                  argv[0], optarg, rundata->nprocs);
                          fprintf(stderr,"  Using %d\n", rundata->npx);
                      } else {
                          rundata->npx = tempint;
                          rundata->npy = rundata->nprocs / tempint;
                      }
                      break;

            case 'Y': tempint = atoi(optarg);
                      if (tempint < 1 || tempint > rundata->nprocs) {
                          fprintf(stderr,
                                  "%s: Cannot use number of processors in y direction %s;\n",
                                  argv[0], optarg);
                          fprintf(stderr,"  Using %d\n", rundata->npy);
                      } else if (rundata->nprocs % tempint != 0) {
                          fprintf(stderr,
                                  "%s: Number of processors in y direction %s does not divide %d;\n",
                                  argv[0], optarg, rundata->nprocs);
                          fprintf(stderr,"  Using %d\n", rundata->npy);
                      } else {
                          rundata->npy = tempint;
                          rundata->npx = rundata->nprocs / tempint;
                      }
                      break;

            case 'f': strncpy(rundata->filename, optarg, MAXFILENAME-1);
                      if (!(tst=fopen(rundata->filename,"w"))) {
                          fprintf(stderr,
                                  "%s: Cannot use filename %s;\n",
                                  argv[0],rundata->filename);
                          fprintf(stderr, "  Using %s\n",defaultfname);
                          strcpy(rundata->filename, defaultfname);
                      } else
                          fclose(tst);
                      break; 

            case 'h':
                  puts("Options: ");
                  puts("    --nx=N         (-x N): Set the number of grid cells in x direction.");
                  puts("    --ny=N         (-y N): Set the number of grid cells in y direction.");
                  puts("    --npx=N        (-X N): Set the number of processors in the x direction.");
                  puts("    --npy=N        (-Y N): Set the number of processors in the y direction.");
                  puts("    --filename=S   (-f S): Set the output filename.");
                  puts("");
                  return +1;

            default: printf("Invalid option %s\n", optarg);
                 break;
        }
    }

    rundata -> myy = rundata->rank / (rundata->npx);
    rundata -> myx = rundata->rank % (rundata->npx);
 
    rundata->localnx = (rundata->globalnx / rundata->npx);
    rundata->localny = (rundata->globalny / rundata->npy);

    /* last row/column gets any extra / fewer points to make things work out: */
    if (rundata->myx == rundata->npx-1) 
        rundata->localnx = (rundata->globalnx - (rundata->npx-1)*(rundata->localnx));
    if (rundata->myy == rundata->npy-1) 
        rundata->localny = (rundata->globalny - (rundata->npy-1)*(rundata->localny));

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

void fillarray2d(double **d, rundata_t *r) {
    int i,j;
    double r2,r2x;
    int gnx = r->globalnx;
    int gny = r->globalny;
    int nx = r->localnx;
    int ny = r->localny;
    int npx = r->npx;
    int npy = r->npy;
    int myx = r->myx;
    int myy = r->myy;
    int startx,starty;
    double sigma=gnx/4.;

    startx = (gnx/npx)*myx;
    starty = (gny/npy)*myy;

    for (i=0;i<nx;i++) {
        r2x = ((i+startx)-gnx/2.)*((i+startx)-gnx/2.);
        for (j=0;j<ny;j++) {
            r2 = r2x + ((j+starty)-gny/2.)*((j+starty)-gny/2.);
            d[i][j] = 1. + 4.*exp(-r2/(2.*sigma*sigma));
        }
    } 

    return;
}

void fillarray3d(double ***v, rundata_t *r) {
    int i,j;
    int gnx = r->globalnx;
    int gny = r->globalny;
    int nx = r->localnx;
    int ny = r->localny;
    int npx = r->npx;
    int npy = r->npy;
    int myx = r->myx;
    int myy = r->myy;
    int startx,starty;
    double sigma=gnx/4.;
    double r2,r2x;

    startx = (gnx/npx)*myx;
    starty = (gny/npy)*myy;

    for (i=0;i<nx;i++) {
        r2x = ((i+startx)-gnx/2.)*((i+startx)-gnx/2.);
        for (j=0;j<ny;j++) {
            r2 = r2x + ((j+starty)-gny/2.)*((j+starty)-gny/2.);
            if (r2 == 0) {
                v[0][i][j] = 0.;
                v[1][i][j] = 0.;
            } else {
                v[0][i][j] = exp(-r2/(2.*sigma*sigma))*((j+starty)-gny/2.);
                v[1][i][j] = -exp(-r2/(2.*sigma*sigma))*((i+startx)-gnx/2.);
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

void nearsquare(int nprocs, int *npx, int *npy) {
    int sq = ceil(sqrt((double)nprocs));
    int n,m;

    if (sq*sq == nprocs) {
        *npx = *npy = sq;
    } else {
        for (n = sq+1; n>=1; n--) {
            if (nprocs % n == 0) {
                m = nprocs/n;
                *npx = n; *npy = m;
                if (m<n) {*npx = m; *npy = n;}
                break;
            } 
        }
    }
    return;
}
