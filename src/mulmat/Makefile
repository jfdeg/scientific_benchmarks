# Makefile for gesvd example

# Compiler
CC=gcc
NVCC=nvcc
CFLAGS=-O3 -march=native

# Intel MKL library
MKL_LIB=${MKL_DIR}/lib/intel64
MKL_INC=${MKL_DIR}/include

OPENBLAS_LIB=${OPENBLAS_DIR}/lib/
OPENBLAS_INC=${OPENBLAS_DIR}/include/

MKL_SEQ=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
MKL_OMP=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread 

# Nvidia CUDA library
CARMA_LIB=/usr/local/cudarm4/cuda/lib
CUDA_LIB = -lcublas -lcudart 
CUDA_INCLUDE =/usr/local/cuda/include
# you have to specify the compute capability of your GPU
CUDA_ARCH = -arch=sm_61

# Other libraries
LIB_STD=-lrt -lm

usage:
	@echo "Command available:"
	@echo ""
	@echo "    build double precision (fp64) GPU benchmark :"
	@echo "    make cuda_gemm_fp64"
	@echo "    build single precision (fp32) GPU bencharmk :"
	@echo "    make cuda_gemm_fp32"
	@echo "    build double precision (fp64) CPU bencharmk singlethread) :"
	@echo "    make mkl_gemm_fp64_monothread"
	@echo "    build single precision (fp32) CPU bencharmk singlethread) :"
	@echo "    make mkl_gemm_fp32_monothread"
	@echo "    make openblas_gemm_fp32_monothread"
	@echo "    build double precision (fp64) CPU bencharmk multithreads) :"
	@echo "    make mkl_gemm_fp64_multithread"
	@echo "    build single precision (fp32) CPU bencharmk multithreads) :"
	@echo "    make mkl_gemm_fp32_multithread"
	@echo ""
	@echo "    make clean all builds"
	@echo "Feel free to send me all your results ;)"

cuda_gemm_fp64:
	${NVCC} -m64 $(CUDA_ARCH) -O3  -o ./bin/cuda_gemm_fp64.x ./src/zgemm.cu -I${CUDA_INCLUDE} $(CUDA_LIB) $(LIB_STD)

cuda_gemm_fp32:
	${NVCC} -m64 $(CUDA_ARCH) -O3 -o ./bin/cuda_gemm_fp32.x ./src/cgemm.cu -I${CUDA_INCLUDE} $(CUDA_LIB) $(LIB_STD)

mkl_gemm_fp64_monothread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/mkl_gemm_fp64_monothread.x ./src/zgemm_MKL.c -I$(MKL_INC) -L${MKL_LIB}  $(MKL_SEQ) $(LIB_STD)

mkl_gemm_fp64_multithread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/mkl_gemm_fp64_multithread.x ./src/zgemm_MKL.c -I$(MKL_INC) -L${MKL_LIB}  $(MKL_OMP) $(LIB_STD)

mkl_gemm_fp32_multithread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/mkl_gemm_fp32_multithread.x ./src/cgemm_MKL.c -I$(MKL_INC)  -L${MKL_LIB}  $(MKL_OMP) $(LIB_STD)

mkl_gemm_fp32_monothread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/mkl_gemm_fp32_monothread.x ./src/cgemm_MKL.c -I$(MKL_INC)  -L${MKL_LIB}  $(MKL_SEQ) $(LIB_STD)
	
atlas_gemm_fp32_monothread:
	${CC} -fopenmp $(CFLAGS) -o ./bin/atlas_gemm_fp32_monothread.x ./src/cgemm_ATLAS.c -I${ATLAS_INC} -L${ATLAS_LIB} -lcblas -latlas  $(LIB_STD)	

openblas_gemm_fp32:
	${CC} -fopenmp $(CFLAGS) -o ./bin/openblas_gemm_fp32.x ./src/cgemm_OPENBLAS.c -I${OPENBLAS_INC} -L${OPENBLAS_LIB} -lopenblas  -lpthread $(LIB_STD)

openblas_gemm_fp64:
	${CC} -fopenmp $(CFLAGS) -o ./bin/openblas_gemm_fp64.x ./src/zgemm_OPENBLAS.c -I${OPENBLAS_INC} -L${OPENBLAS_LIB} -lopenblas  -lpthread $(LIB_STD)

clean:
	rm -f *.x *.exe bin/*.x

