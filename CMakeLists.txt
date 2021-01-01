cmake_minimum_required(VERSION 3.11)

project(scientific_benchmark LANGUAGES C CXX)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/.cmake/")

# find_package(CUDA)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(GTest REQUIRED)

enable_testing()
include(GoogleTest)


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

set(COVMAT_MKL_SRC src/covmat/mkl/covmat_mkl_benchmark_sp.c src/covmat/mkl/covmat_mkl_complex_sp.c include/covmat_mkl.h)
set(MKL_LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -lrt -lpthread -ldl)
set(C_FLAGS -g -fopenmp -march=native -fno-stack-protector -flto -funroll-loops -ffast-math -Wall)

file(GLOB TEST_SRC_FILES tests/src/*.cpp)
file(GLOB INT_TEST_SRC_FILES tests/src/integration/*.cpp)
file(GLOB UTILS_SRC_FILES src/utils/*.c)

add_executable(covmat_mkl_benchmark_sp src/covmat/main.c ${COVMAT_MKL_SRC} ${UTILS_SRC_FILES})
target_link_libraries(covmat_mkl_benchmark_sp ${MKL_LIBS})
target_include_directories(covmat_mkl_benchmark_sp PRIVATE include)
target_compile_options(covmat_mkl_benchmark_sp PUBLIC ${C_FLAGS})

add_executable(tests ${TEST_SRC_FILES} ${UTILS_SRC_FILES} ${COVMAT_MKL_SRC})
target_include_directories(tests PRIVATE include)
add_dependencies(tests covmat_mkl_benchmark_sp)
target_link_libraries(tests GTest::GTest ${MKL_LIBS})
gtest_discover_tests(tests)

add_executable(integration_tests ${INT_TEST_SRC_FILES} ${TEST_SRC_FILES} ${UTILS_SRC_FILES} ${COVMAT_MKL_SRC})
target_include_directories(integration_tests PRIVATE include)
add_dependencies(integration_tests covmat_mkl_benchmark_sp)
target_link_libraries(integration_tests GTest::GTest ${MKL_LIBS})
gtest_discover_tests(integration_tests)