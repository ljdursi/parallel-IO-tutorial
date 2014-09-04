#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define FILESIZE 1024*1024*128

int write_file_bin(const char *fname, const double *data, const int ndata) {

    FILE *fp;
    time_t start, end;

    fp=fopen(fname,"wb");
    assert(fp);
    start = time(NULL);
    fwrite(data, sizeof(double), ndata, fp);
    end = time(NULL);
    fclose(fp);

    return (int)(end-start);
}

int write_file_ascii(const char *fname, const double *data, const int ndata) {

    FILE *fp;
    time_t start, end;
    int i;

    fp=fopen(fname,"wb");
    assert(fp);
    start = time(NULL);
    for (i=0;i<ndata;i++) {
        fprintf(fp,"%lf\n",data[i]);
    }
    end = time(NULL);
    fclose(fp);

    return (int)(end-start);
}

int main(int argc, char **argv) {
    double *data;
    int i;
    int asciitime, bintime;

    data = (double *)malloc(FILESIZE * sizeof(double));
    assert(data);
    for (i=0;i<FILESIZE;i++) {
        data[i] = i*(double)i/2.;
    }

    asciitime = write_file_ascii("data.txt",data,FILESIZE); 
    bintime   = write_file_bin("data.dat",data,FILESIZE); 

    printf("Time to write files: ASCII: %d, Binary: %d\n",asciitime, bintime);

    return 0;
}

