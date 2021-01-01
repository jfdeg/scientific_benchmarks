#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include "utils.h"
#include "benchmarks.h"
#include "covmat_mkl.h"

benchmark_results covmat_mkl_benchmark_sp(int matrix_size, int n_runs) {
	MKL_Complex8 *signal_in  = (MKL_Complex8*) malloc(matrix_size * matrix_size * sizeof(MKL_Complex8));
	MKL_Complex8 *signal_out = (MKL_Complex8*) malloc(matrix_size * matrix_size * sizeof(MKL_Complex8));

	if (signal_in == NULL || signal_out == NULL) {
		fprintf(stderr, "[Error]: memory allocation failed, occured in %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	benchmark_results results = covmat_mkl_complex_sp(signal_out, signal_in, matrix_size, matrix_size, n_runs);

	free(signal_in);
	free(signal_out);

	return results;
}
