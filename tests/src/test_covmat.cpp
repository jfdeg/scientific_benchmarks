#include "gtest/gtest.h"

extern "C" {
#include "benchmarks.h"
#include "covmat_mkl.h"
}

TEST(CovMatrixMKL_SP, basic) {
    int matrix_size = 100, n_runs = 1;
    benchmark_results r = covmat_mkl_benchmark_sp(matrix_size, n_runs);
}

TEST(CovMatrixMKL_SP, 10Runs) {
    int matrix_size = 100, n_runs = 10;
    benchmark_results r = covmat_mkl_benchmark_sp(matrix_size, n_runs);
}