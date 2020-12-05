## scientific_benchmarks
* [General info](#general-info)
* [Setup](#setup)

## General info
This repo contains several benchmarks that can measure CPU and GPU performance. 

CPU code either uses Intel MKL lib (optimized for Intel x86 CPU but also runs on x86 AMD CPU) or OpenBLAS (can ran on all architectures). 

GPU code either uses CUDA libs (only works with Nvidia GPUs) or OpenCL

List of available benchmarks :

* CovarianceMatrixEstimation : computes the covariance matrix of N vectors of size N using a matrix-matrix multiplication (gemm)
	
	
## Setup
First, you have to generate the data required by the benchmarks by using the data_generation MATLAB script located in the data folder (Python version will come later). 

To run the benchmarks, please read the README included in each foler.



