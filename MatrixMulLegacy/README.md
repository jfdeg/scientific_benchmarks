# gemm_benchmark

Simple program that computes complex matrix-matrix multiplication on CPU et and on GPU in both FP32/FP64 mode.

CPU benchmark uses Intel Math Kernel Library and OpenBLAS

GPU benchmark uses CUDA

You should edit the Makefile to compile binaries


Binaries environment variables:
MATMUL_M : constant M value
MATMUL_N : constant N value
MATMUL_K : constant K value
TRANSA   : Transpose matrix A. choose between N or T (default N)
TRANSB   : Transpose matrix B. choose between N or T (default N)