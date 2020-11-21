#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */
#include <stdio.h>
#include <stdlib.h>

//masquer les PRINTF en mode release
//addcsv donne une sortie formatee pour les csv
#ifndef RELEASE
#define PRINTF printf
#define ADDCSV
#else
#define PRINTF
#define ADDCSV printf
#endif

extern void benchCovMatEstimation_MKL_dp();
extern void benchCovMatEstimation_MKL_sp();

#ifdef CL
extern void benchCovMatEstimation_OPENCL_sp();
#endif

#ifdef CUDA
extern void benchCovMatEstimation_CUDA_dp();
extern void benchCovMatEstimation_CUDA_sp();
#endif

int main(int argc, char **argv) {



	PRINTF("\n\n");
	PRINTF("====================\t CovMatEstimation (x%d)\t ====================\n\n", nb_fois);
#ifdef CL
	PRINTF("Starting OPENCL Simple Precision...\n\n");
	ADDCSV("CovMatEstimation;CL_SP;FO;%d;",type_fo);
	benchCovMatEstimation_OPENCL_sp();
#endif
	PRINTF("Starting MKL Double Precision...\n\n");
	ADDCSV("CovMatEstimation;MKL_DP;FO;%d;",type_fo);
	benchCovMatEstimation_MKL_dp();
	PRINTF("\n\n");

	PRINTF("Starting MKL Simple Precision...\n\n");
	ADDCSV("CovMatEstimation;MKL_SP;FO;%d;",type_fo);
	benchCovMatEstimation_MKL_sp();
	PRINTF("\n\n");
#ifdef CUDA
	PRINTF("Starting CUDA Double Precision...\n\n");
	ADDCSV("CovMatEstimation;CUDA_DP;FO;%d;",type_fo);
	benchCovMatEstimation_CUDA_dp();
	PRINTF("\n\n");

	PRINTF("Starting CUDA Simple Precision...\n\n");
	ADDCSV("CovMatEstimation;CUDA_SP;FO;%d;",type_fo);
	benchCovMatEstimation_CUDA_sp();
	PRINTF("\n\n");

#endif
	return EXIT_SUCCESS;
}

