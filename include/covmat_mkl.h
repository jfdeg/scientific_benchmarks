#pragma once

#include "benchmarks.h"
#ifndef OPEN
#include <mkl.h>
#endif

benchmark_results covmat_mkl_benchmark_sp(int matrix_size, int n_runs);
benchmark_results covmat_mkl_complex_sp(MKL_Complex8 *signal_in, MKL_Complex8 *signal_out, int N, int M, int n_runs);
