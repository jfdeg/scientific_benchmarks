
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include <mkl.h>


void readData(MKL_Complex8 *signal, size_t size_sig, MKL_Complex8 *reference, size_t size_ref) {
	MKL_Complex16 *tmp = (MKL_Complex16*) malloc(size_ref*sizeof(MKL_Complex16));

	if (tmp == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

    // a mettre a jour .....

    char path[150];
    char *home = "/home/jfd/projets/scientific_benchmarks/data/CovMatEstimation/";

	strcpy(path, home);
	strcat(path, "CME_Complex_out_100.bin");

	FILE *fid = fopen(path, "rb");
	if (fid == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	size_t count = fread(tmp, sizeof(MKL_Complex16), size_ref, fid);

	if (count < (size_t) size_ref) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	for(size_t i = 0 ; i < (size_t) size_ref ; i++) {
		(*(reference+i)).real = (float) (*(tmp+i)).real;
		(*(reference+i)).imag = (float) (*(tmp+i)).imag;
	}

	free(tmp);
	fclose(fid);

	tmp = (MKL_Complex16*) malloc(size_sig * sizeof(MKL_Complex16));
	if (tmp == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	strcpy(path, home);
	strcat(path, "CME_Complex_in_100.bin");

	fid = fopen(path, "rb");
	if (fid == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	count = fread(tmp , sizeof(MKL_Complex16), size_sig, fid);
	if (count < (size_t) size_sig) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	for (size_t i = 0 ; i < (size_t) size_sig ; i++) {
		(*(signal+i)).real = (float) (*(tmp+i)).real;
		(*(signal+i)).imag = (float) (*(tmp+i)).imag;
	}

	fclose(fid);
}