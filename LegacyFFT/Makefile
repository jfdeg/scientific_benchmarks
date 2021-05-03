
# Compiler
CC=gcc
NVCC=nvcc
CFLAGS=-O3 -march=native

# Intel MKL library
MKL_LIB=${MKL_DIR}/lib/intel64
LIB_OMP=${MKL_ROOT}/compilers_and_libraries/linux/lib/intel64
MKL_INC=${MKL_DIR}/include

MKL_SEQ=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
MKL_OMP=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Nvidia CUDA library
CARMA_LIB=/usr/local/cudarm4/cuda/lib
CUDA_LIB = -lcublas -lcudart -lcufft
CUDA_INCLUDE =/usr/local/cuda/include
# you have to specify the compute capability of your GPU
CUDA_ARCH = -arch=sm_61

# Other libraries
LIB_STD=-lrt -lm

usage:
	@echo "Command available:"
	@echo ""
	@echo "    build double precision (fp64) GPU benchmark :"
	@echo "    make cuda_fft_fp64"
	@echo "    build single precision (fp32) GPU bencharmk :"
	@echo "    make cuda_fft_fp32"
	@echo "    build single precision (fp32) CPU bencharmk singlethread) :"
	@echo "    make mkl_fft_fp32_monothread"
	@echo "    build double precision (fp64) CPU bencharmk multithreads) :"
	@echo "    make mkl_fft_fp32_multithread"
	@echo ""
	@echo "    make clean all builds"
	@echo "Feel free to send me all your results ;)"

cuda_fft_fp64:
	${NVCC} -m64 $(CUDA_ARCH)  -O3  -o ./bin/cuda_fft_fp64.x ./src/zgemm.cu -I${CUDA_INCLUDE} $(CUDA_LIB) $(LIB_STD)

cuda_fft_fp32:
	${NVCC} -m64 $(CUDA_ARCH)  -O3 -o ./bin/cuda_fft_fp32.x ./src/cuFFT_32_flops.cu -I${CUDA_INCLUDE} $(CUDA_LIB) $(LIB_STD)

mkl_fft_fp32_multithread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/mkl_fft_fp32_multithread.x ./src/MKL_flops_32_omp.c -I$(MKL_INC)  -L${MKL_LIB}  $(MKL_OMP) $(LIB_STD)

mkl_fft_fp32_monothread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/mkl_fft_fp32_monothread.x ./src/MKL_flops_32_omp.c -D MKL_LIB -DNBCORES=1 -I$(MKL_INC) -L${MKL_LIB}  $(MKL_SEQ) $(LIB_STD)
	

clean:
	rm -f *.x *.exe

