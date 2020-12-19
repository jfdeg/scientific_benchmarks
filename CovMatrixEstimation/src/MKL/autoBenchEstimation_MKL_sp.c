#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#ifndef RELEASE
#define PRINTF printf
#define ADDCSV
#else
#define PRINTF
#define ADDCSV printf
#endif

// Includes spécifique à la librairie MKL utilisée
#ifndef OPEN
#include <mkl.h>
#endif

#ifndef USE_MKL_MALLOC
	#define MALLOC(a) malloc(a)
	#define FREE(a) free(a)
#else
	#ifdef OPEN
		#define MALLOC(a) malloc(a) //aligned_alloc(a,64)
		#define FREE(a)   free(a)
	#else
		#define MALLOC(a) mkl_malloc(a,64)
		#define FREE(a) mkl_free(a)
	#endif
#endif


extern int mklEstimation_Complex_sp(MKL_Complex8 *signal_out, MKL_Complex8 *signal_in, int N, int M, int nbMat, char transN, char transC, MKL_Complex8 norm, MKL_Complex8 null);
extern void readData(MKL_Complex8 *signal, size_t size_sig, MKL_Complex8 *reference, size_t size_ref);

void AutoBenchEstimation_MKL_sp(int nbRUN) {

    /* declarations */
    size_t size_mat, size_sig, size_ref;
    float total,min_time,max_time,elapsedTime,nb_operation;
	struct timespec tpdeb;
	struct timespec tpfin;
	clockid_t clock_id = CLOCK_REALTIME;
    char transN, transC;
    MKL_Complex8 norm, null, diff;
	float error, sum;
    

    /* Loop on matrix size */
    for(int N = 100 ; N <= 100 ; N = N + 100) {

		/** Data init **/
		size_mat = (size_t) (N * N);
		size_sig = (size_t) (size_mat * nbRUN);
		size_ref = (size_t) (size_mat * nbRUN);

		
		MKL_Complex8 *signal_in = (MKL_Complex8*) MALLOC(size_sig*sizeof(MKL_Complex8));
		MKL_Complex8 *signal_out = (MKL_Complex8*) MALLOC(size_ref*sizeof(MKL_Complex8));

		if (signal_in == NULL || signal_out == NULL) {
			fprintf(stderr, "[Error]: memory allocation failed, occured in %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
		}

		/* Timing tools  ********************************************************************************************/
		total = 0.0f;
		min_time = 65536.0f;
		max_time = 0.0f;
		elapsedTime  = 0.0f;
		nb_operation = (float) (1e-9 * N * N * 8.0f);

		/* gemm parameters */
		transN = 'N';
		transC = 'C';
		norm.real = 1.0f/(float) N;
		norm.imag = 0.0f;
		null.real = 0.0f;
		null.imag = 0.0f;

		// Warmup
		memset(signal_out, 0, size_ref*sizeof(MKL_Complex8));
		mklEstimation_Complex_sp(signal_out, signal_in, N, N, 1, transN, transC, norm, null);
		memset(signal_out, 0, size_ref*sizeof(MKL_Complex8));

		// compute and mesure performance
		for(int i = 0 ; i < nbRUN ; i++) {
			clock_gettime(clock_id, &tpdeb);

			mklEstimation_Complex_sp(signal_out + i * size_mat, signal_in + i * size_mat, N, N, 1, transN, transC, norm, null);

			clock_gettime(clock_id, &tpfin);

			elapsedTime = (float) (tpfin.tv_sec - tpdeb.tv_sec) + (float) ((float) (tpfin.tv_nsec - tpdeb.tv_nsec) * 1e-9);

			total += elapsedTime;

			max_time = fmaxf(max_time, elapsedTime);
			min_time = fminf(min_time, elapsedTime);
		}

		PRINTF("*************\tTotal MKL Simple Precision\t*************\n");
		PRINTF("Average elapsed time\t= %lf ms\n", (total/nbRUN)*1000);
		PRINTF("Min time \t= %lf ms\n", min_time*1000);
		PRINTF("Max time \t= %lf ms\n", max_time*1000);
		PRINTF("Nb operation\t= %g Gflop\n", nb_operation);
		PRINTF("Performance\t= %g Gflop/s\n", nb_operation/(total/nbRUN));


		/* free memory **********************************************************************************************/
		FREE(signal_in);
		FREE(signal_out);
	}
}
