#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <math.h>

#include <mkl.h>

#include "utils.h"
#include "benchmarks.h"

/**
 * @brief Estimate of the covariance matrix using cgemm
 * 
 * @param signal_in Input matrix of shape (n, m) and type MKL_Complex8,
 * must be pre-allocated and freed by the caller
 * @param signal_out Output matrix of shape (n, m) and type MKL_Complex8,
 * must be pre-allocated and freed by the caller
 * 
 */
benchmark_results covmat_mkl_complex_sp(
	MKL_Complex8 *signal_in, MKL_Complex8 *signal_out, int n, int m, int n_runs
) {
    /* time vars and consts */
	struct timespec tpdeb;
	struct timespec tpfin;
	static const clockid_t clock_id = CLOCK_REALTIME;
    float duration = 0.0f;
    float min_time = 65536.0f;
    float max_time = 0.0f;
    float total_time  = 0.0f;
    int matrix_size = n * m;
    float gflop = (float) (1e-9 * matrix_size * 8.0f);

    /* gemm parameters */
    static const char transN = 'N', transC = 'C';
    MKL_Complex8 norm = {1.0f/(float) matrix_size, norm.imag = 0.0f};
    MKL_Complex8 beta = {0.0f, 0.0f};

    // Warmup
    memset(signal_out, 0, matrix_size * sizeof(MKL_Complex8));
    cgemm(&transN, &transC, &n, &n, &m, &norm, signal_in, &n, signal_in, &n, &beta, signal_out, &n);

    // compute and mesure performance
    for (int i = 0 ; i < n_runs ; i++) {
        clock_gettime(clock_id, &tpdeb);

        cgemm(&transN, &transC, &n, &n, &m, &norm, signal_in, &n, signal_in, &n, &beta, signal_out, &n);

        clock_gettime(clock_id, &tpfin);

        duration = (tpfin.tv_sec - tpdeb.tv_sec) + (tpfin.tv_nsec - tpdeb.tv_nsec) * 1e-9;

        total_time += duration;

        max_time = fmaxf(max_time, total_time);
        min_time = fminf(min_time, total_time);
    }

	float avg_time = total_time / n_runs;
	float gflops = gflop / avg_time;

    // TODO: find a better way to format extra_infos and load it in the struct
    char extra_infos[128];
    sprintf(extra_infos, "matrix size: %d x %d, n_runs: %d", n, m, n_runs);

	benchmark_results bench_results = {
		"Covariance matrix estimation MKL",
        "",
		n_runs,
		total_time,
		min_time,
		avg_time,
		max_time,
		gflop,
		gflops,
        SP
	};

    memcpy(&bench_results.extra_infos, extra_infos, 128);

	return bench_results;
}