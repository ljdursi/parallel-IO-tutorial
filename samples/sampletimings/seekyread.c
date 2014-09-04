#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define FILESIZE 1024*1024*128

int read_bin_onego(const char *fname, double *data, const int ndata) {

    FILE *fp;
    time_t start, end;

    fp=fopen(fname,"rb");
    assert(fp);
    start = time(NULL);
    fread(data, sizeof(double), ndata, fp);
    end = time(NULL);
    fclose(fp);

    return (int)(end-start);
}

int read_bin_chunky(const char *fname, double *data, const int ndata) {

    FILE *fp;
    time_t start, end;
    const int CHUNKSIZE=1024;
    int nchunks = FILESIZE/CHUNKSIZE;

    fp=fopen(fname,"rb");
    assert(fp);
    start = time(NULL);
    for (int i=0; i<nchunks; i++) {
        int start = i*CHUNKSIZE;
        fread(&(data[start]), sizeof(double), CHUNKSIZE, fp);
    }
    end = time(NULL);
    fclose(fp);

    return (int)(end-start);
}

int read_bin_seeky(const char *fname, double *data, const int ndata) {

    FILE *fp;
    time_t start, end;
    const int CHUNKSIZE=1024;
    int nchunks = FILESIZE/CHUNKSIZE;

    fp=fopen(fname,"rb");
    assert(fp);
    start = time(NULL);
    for (int i=0; i<nchunks; i++) {
        int chunk = (i*13 % nchunks);
        int start = chunk*CHUNKSIZE;
        fseek(fp, start*sizeof(double), SEEK_SET);
        fread(&(data[start]), sizeof(double), CHUNKSIZE, fp);
    }
    end = time(NULL);
    fclose(fp);

    return (int)(end-start);
}

int read_bin_iopsy(const char *fname, double *data, const int ndata) {

    FILE *fp;
    time_t start, end;
    const int CHUNKSIZE=1024;
    int nchunks = FILESIZE/CHUNKSIZE;

    assert(fp);
    start = time(NULL);
    for (int i=0; i<nchunks; i++) {
        int chunk = (i*13 % nchunks);
        int start = chunk*CHUNKSIZE;
        fp=fopen(fname,"rb");
        fseek(fp, start*sizeof(double), SEEK_SET);
        fread(&(data[start]), sizeof(double), CHUNKSIZE, fp);
        fclose(fp);
    }
    end = time(NULL);

    return (int)(end-start);
}

int main(int argc, char **argv) {
    double *data;
    int i;
    int contigtime, chunktime, seektime, openclose;

    data = (double *)malloc(FILESIZE * sizeof(double));
    assert(data);

    chunktime  = read_bin_chunky("data.dat",data,FILESIZE); 
    contigtime = read_bin_onego("data.dat",data,FILESIZE); 
    seektime   = read_bin_seeky("data.dat",data,FILESIZE); 
    openclose  = read_bin_iopsy("data.dat",data,FILESIZE); 

    printf("Time to read files: One go: %d, Chunky: %d, Seeky: %d, iopsy: %d\n",contigtime, chunktime, seektime, openclose);

    return 0;
}

