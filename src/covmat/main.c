#include <stdio.h>
#include <stdlib.h>

#include "covmat_mkl.h"
#include "io.h"


int main(int argc, char **argv) {

	int matrix_size = 1000;
    int n_runs = 5;
	benchmark_results r;

	puts("Starting benchmarks...\n");

	r = covmat_mkl_benchmark_sp(matrix_size, n_runs);
	display_benchmark_results(r);

	puts("");

	return EXIT_SUCCESS;
}
