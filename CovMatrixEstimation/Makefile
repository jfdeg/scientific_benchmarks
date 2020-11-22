# Locati################################################################################################################################
# Nom........... : Makefile
# Role.......... : Compile le MatrixMul benchmark 
# Auteur........ : J-F DEGURSE, F. TURBAN, A. DEGURSE 
# Date.......... : 20/11/2020
################################################################################################################################

# Compilateur C
CC   ?= gcc
CL   ?= 0
CUDA ?= 0
OPEN ?= 0

# Architecture NVIDIA
NV_ARCH=-gencode arch=compute_70,code=sm_70 # TITAN V => 7.0

# Listes des fichiers .c et .o pour MKL double precision
#SRC_C_DP = $(wildcard src/MKL/*_dp.c)
#OBJ_C_DP = $(subst src, obj, $(SRC_C_DP:.c=.o))

# Listes des fichiers .c et .o pour MKL simple precision
SRC_C_SP = $(wildcard src/MKL/*_sp.c )
OBJ_C_SP = $(subst src, obj, $(SRC_C_SP:.c=.o))

# Listes des fichiers .c et .o pour OpenCL simple precision
ifeq ($(CL),1)
SRC_C_SP += $(wildcard src/CL/*_sp.c )
OBJ_C_SP = $(subst src, obj, $(SRC_C_SP:.c=.o))
CFLAGS += -DCL
#OBJ_C_SP = $(subst obj/CL, obj,$(OBJ_C_SP_tmp))
endif 

# Listes des fichiers .cu et .o pour CUDA simple precision
ifeq ($(CUDA),1)
SRC_CUDA_SP = $(wildcard src/CUDA/*_sp.cu)
OBJ_CUDA_SP = $(subst src, obj, $(SRC_CUDA_SP:.cu=.o))

# Listes des fichiers .cu et .o pour CUDA double precision
SRC_CUDA_DP = $(wildcard src/CUDA/*_dp.cu)
OBJ_CUDA_DP = $(subst src, obj, $(SRC_CUDA_DP:.cu=.o))

CUDA_BENCH_OBJ=obj/CUDA/benchCUDA.o
CFLAGS += -DCUDA
endif

# Flags C 
CFLAGS +=  -std=gnu99 -O3 -Wwrite-strings -Wstrict-prototypes 		\
         -Wuninitialized -Wno-missing-braces -Wno-missing-field-initializers -DUSE_MKL_MALLOC

# Variables associees a MKL
MKL_INCLUDE  = -I$(MKL_INCLUDE_PATH) 

ifeq ($(CC),icc)
MKL_FLAGS    = -openmp -xAVX -axCORE-AVX2 -ipo
MKL_LIBRARY  = -L$(MKL_LIBRARY_PATH) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lm -lrt -lpthread -ldl
else # gcc
MKL_FLAGS    ?=  -fopenmp -march=native -flto -funroll-loops -ffast-math
MKL_LIBRARY  = -L$(MKL_LIBRARY_PATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lm -lrt -lpthread -ldl
endif

# Variables associees a OpenCL
ifeq ($(CL),1)
CL_PATH ?= /data/enf/OpenCL
CL_INCLUDE = -I$(CL_PATH)/include/ 
CL_LIBRARY = -L$(CL_PATH)/lib/ -lclblast -lOpenCL
endif

# Variables associees a CUDA
ifeq ($(CUDA),1)
CUDA_INCLUDE = -I$(CUDA_INCLUDE_PATH) -I$(CUDA_PATH)/samples/common/inc
CUDA_LIBRARY = -L$(CUDA_LIBRARY_PATH) -lcusolver -lcublas -lcufft -lcudart -lrt -lpthread -ldl 
CUDA_FLAGS   += -Wall -Wextra -O3 -fopenmp
endif

ifeq ($(OPEN),1)
OPENBLAS_PATH ?= /data/enf/STAP/openblas/single/
MKL_INCLUDE   = -I$(OPENBLAS_PATH)/include
MKL_LIBRARY   = -L$(OPENBLAS_PATH)/lib/ -Wl,--start-group -lopenblas -Wl,--end-group -lm -lrt -lpthread -ldl
CFLAGS += -DOPEN
endif

# Variables associées à MEX
MEX_FLAGS = -fPIC -shared

all : benchCovMatEstimation_MKL

run : benchCovMatEstimation_MKL
	./benchCovMatEstimation_MKL.out


benchCovMatEstimation : $(OBJ_C_SP) $(OBJ_C_DP) $(OBJ_CUDA_SP) $(OBJ_CUDA_DP) obj/main.o
ifeq ($(CUDA),1)
	nvcc -m64 $(NV_ARCH) --compiler-options "\$(CUDA_FLAGS) \$(MEX_FLAGS)" --device-link $(OBJ_CUDA_SP) $(OBJ_CUDA_DP) -o $(CUDA_BENCH_OBJ)
endif
	$(CC) -o $@.out $+ $(CUDA_BENCH_OBJ) $(MKL_FLAGS) $(MKL_LIBRARY) $(CUDA_FLAGS) $(CUDA_LIBRARY) $(CL_LIBRARY)
	

benchCovMatEstimation_MKL : $(OBJ_C_SP) $(OBJ_C_DP) obj/main.o
	$(CC) -o $@.out $+ $(MKL_FLAGS) $(MKL_LIBRARY) 

ifeq ($(CUDA),1)
obj/%_dp.o : src/%_dp.cu
	nvcc -ccbin=cuda-gcc -m64 $(NV_ARCH) -use_fast_math --compiler-options "\$(CUDA_FLAGS) \$(MEX_FLAGS)" $(CUDA_INCLUDE) --device-c $< -o $@

obj/%_sp.o : src/%_sp.cu
	nvcc -ccbin=cuda-gcc -m64 $(NV_ARCH) -use_fast_math --compiler-options "\$(CUDA_FLAGS) \$(MEX_FLAGS)" $(CUDA_INCLUDE) --device-c $< -o $@
endif 

#obj/%_dp.o : src/%_dp.c
#	$(CC) -o $@ -c $^ $(CFLAGS) $(MKL_FLAGS) $(MKL_INCLUDE)

obj/%_sp.o : src/%_sp.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(MKL_FLAGS) $(MKL_INCLUDE)
	
ifeq ($(CL),1)
obj/CL/%_sp.o : src/CL/%_sp.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(MKL_FLAGS)  $(CL_INCLUDE)	
endif

obj/main.o : src/main.c
	$(CC) $(CFLAGS) $(MKL_FLAGS) $(MKL_INCLUDE) $(CUDA_FLAGS) $(CUDA_INCLUDE) -c $< -o $@
	
clean :
	rm -rf obj/*/*
	rm -rf obj/main.o
	rm -rf *.out
