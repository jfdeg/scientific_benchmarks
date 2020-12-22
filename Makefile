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


# Listes des fichiers .c et .o pour MKL simple precision
COVMAT_SRC_C_SP = $(wildcard src/covmat/mkl/*_sp.c )
COVMAT_OBJ_C_SP = $(subst src, bin, $(COVMAT_SRC_C_SP:.c=.o))

UTILS_SRC_C = $(wildcard src/utils/*.c )
UTILS_OBJ_C = $(subst src, bin, $(UTILS_SRC_C:.c=.o))

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
MKL_LIBRARY  = -L$(MKL_LIBRARY_PATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -lm -lrt -lpthread -ldl
endif


# Variables associées à MEX
MEX_FLAGS = -fPIC -shared

all : gen_build_dir covmat_mkl_benchmark

run : covmat_mkl_benchmark
	./bin/covmat_mkl_benchmark.out

# project include
COVMAT_INCLUDE = -I include

gen_build_dir:
	mkdir -p bin/covmat/mkl bin/utils

covmat_mkl_benchmark : $(COVMAT_OBJ_C_SP) $(UTILS_OBJ_C) bin/covmat/main.o
	mkdir -p bin/covmat/mkl
	$(CC) -o bin/$@.out $+ $(MKL_FLAGS) $(MKL_LIBRARY)

bin/utils/%.o : src/utils/%.c
	$(CC) -o $@ -c $^ $(CFLAGS) -I include

bin/covmat/mkl/%_sp.o : src/covmat/mkl/%_sp.c
	$(CC) -o $@ -c $^ $(CFLAGS) $(MKL_FLAGS) $(MKL_INCLUDE) -I include

bin/covmat/main.o : src/covmat/main.c
	$(CC) $(CFLAGS) $(MKL_FLAGS) $(MKL_INCLUDE) $(CUDA_FLAGS) $(CUDA_INCLUDE) -Iinclude -c $< -o $@
	
clean :
	@rm -rf bin build

cmake:
	cmake -DCMAKE_BUILD_TYPE=Debug -S . -B build/debug && cmake --build build/debug

release:
	cmake -DCMAKE_BUILD_TYPE=release -S . -B build/release && cmake --build build/release

test: cmake
	./build/debug/bin/tests

generate-test-data:
	python3 tests/resources/generate_test_data.py

test-all: cmake
	./build/debug/bin/integration_tests

help: usage

usage:
	@echo "Command available:"
	@echo "    make \t\t- builds the project using this Makefile, test are not supported"
	@echo "    make cmake \t\t- builds the project using this cmake in debug mode, test are built too"
	@echo "    make test \t\t- builds the project using this cmake in debug mode, unit tests are run"
	@echo "    make generate-test-data\t\t- generate test data for integration tests using python3"
	@echo "    make test-all \t- builds the project using this cmake in debug mode, all tests are run"

	@echo "    make clean"