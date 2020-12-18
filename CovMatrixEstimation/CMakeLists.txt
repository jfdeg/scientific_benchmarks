cmake_minimum_required(VERSION 3.0)

project(scientific_benchmark)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/.cmake/")

# find_package(CUDA)

# set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

# if (CUDA_FOUND)
# 	#set(CUDA_HOST_COMPILER /usr/bin/gcc-6)

# 	include_directories(${CUDA_INCLUDE_DIRS})
# 	#link_directories(/opt/cuda/lib64)

# 	set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
# 	set(
# 	    CUDA_NVCC_FLAGS
# 	    ${CUDA_NVCC_FLAGS}
# 	    #-O3 -gencode arch=compute_20,code=sm_20 # became usupported
# 	    -O3 -gencode arch=compute_30,code=sm_30
# 	)


# 	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")


# 	cuda_add_executable(cuda_gemm_fp64 src/zgemm.cu)
# 	target_link_libraries (cuda_gemm_fp64 ${CUDA_CUBLAS_LIBRARIES} )
# 	install(TARGETS cuda_gemm_fp64 DESTINATION bin)


# endif()

# assuming that MKL was installed with a package manager (as in the docker image)
include_directories("/usr/include/mkl")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")

add_executable(autoBenchEstimation_MKL_sp src/MKL/autoBenchEstimation_MKL_sp.c src/MKL/mklEstimation_Complex_sp.c src/main.c)
target_link_libraries(autoBenchEstimation_MKL_sp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -lrt -lpthread -ldl)
