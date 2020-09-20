## scientific_benchmarks
* [General info](#general-info)
* [Setup](#setup)

## General info
This repo contains several benchmarks that can measure CPU and GPU performance. 

CPU code either uses Intel MKL lib (optimized for Intel x86 CPU but also runs on x86 AMD CPU) or OpenBLAS (can ran on all architectures). 

GPU code either uses CUDA libs (only works with Nvidia GPUs) or OpenCL

List of available benchmarks :

* MatrixMul : computes a complex matrix-matrix multiplication using the gemm function.
	
	
## Setup
To run the benchmarks, please read the README included in each foler.



