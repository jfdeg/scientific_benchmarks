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
static void readData(MKL_Complex8 *signal, size_t size_sig, MKL_Complex8 *reference, size_t size_ref);

void AutoBenchEstimation_MKL_sp(int nbRUN) {

    /* declarations */
    size_t size_mat,size_sig,size_ref;
    float total,min_time,max_time,elapsedTime,nb_operation;
	struct timespec tpdeb;
	struct timespec tpfin;
	clockid_t clock_id = CLOCK_REALTIME;
    char transN, transC;
    MKL_Complex8 norm,null,diff;
	float error,sum;
    

    /* Loop on matrix size */
    for(int N = 100 ; N <= 100 ; N = N + 100) {

	/** Data init **/
    size_mat = (size_t) (N*N);
	size_sig = (size_t) (size_mat*nbRUN);
	size_ref = (size_t) (size_mat*nbRUN);

	MKL_Complex8 *signal_in = (MKL_Complex8*) MALLOC(size_sig*sizeof(MKL_Complex8));
	if(signal_in == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	MKL_Complex8 *signal_out = (MKL_Complex8*) MALLOC(size_ref*sizeof(MKL_Complex8));
	if(signal_out == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	MKL_Complex8 *reference = (MKL_Complex8*) MALLOC(size_ref*sizeof(MKL_Complex8));
	if(reference == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	readData(signal_in, size_sig, reference, size_ref);

	/*************************************************************************************************************************/

	/* Timing tools  ********************************************************************************************/
	total = 0.0f;
	min_time = 65536.0f;
	max_time = 0.0f;
	elapsedTime  = 0.0f;
	nb_operation = (float) (1e-9*N*N*8.0f);

	FILE *fid = fopen("MatrixMul_MKL_sp.txt", "w+");
	if(fid == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
		return;
	}

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

	for(int i = 0 ; i < nbRUN ; i++) {
		//memset(signal_out, 0, size_ref*sizeof(MKL_Complex8));

		clock_gettime(clock_id, &tpdeb);

		mklEstimation_Complex_sp(signal_out+i*size_mat, signal_in+i*size_mat, N, N, 1, transN, transC, norm, null);

		clock_gettime(clock_id, &tpfin);

		elapsedTime = (float) (tpfin.tv_sec - tpdeb.tv_sec) + (float) ((float) (tpfin.tv_nsec - tpdeb.tv_nsec) * 1e-9);
		total += elapsedTime;
		max_time = fmaxf(max_time, elapsedTime);
		min_time = fminf(min_time, elapsedTime);
		fprintf(fid, "%g\n", elapsedTime*1000);
	}

	fclose(fid);



	/*************************************************************************************************************************/

	/* Results comparison with MATLAB  *****************************************************************************************/

	error = 0.0f;
	sum   = 0.0f;

	for(int i = 0 ; i < nbRUN ; i++) {
		for(int j = 0 ; j < N ; j++) {
			for(int k = 0 ; k < N ; k++) {
				if(j >= k) {
					MKL_Complex8 *ref = (MKL_Complex8*) reference+k+(j+i*N)*N;
					MKL_Complex8 *sig = (MKL_Complex8*) signal_out+k+(j+i*N)*N;
					diff.real  = ref->real - sig->real;
					diff.imag  = ref->imag - sig->imag;
					error += sqrtf(diff.real * diff.real + diff.imag * diff.imag);
					sum   += sqrtf(ref->real * ref->real + ref->imag * ref->imag);
				}
			}
		}
	}

	PRINTF("*************\tTotal MKL Simple Precision\t*************\n");
	PRINTF("Average elapsed time\t= %lf ms\n", (total/nbRUN)*1000);
	PRINTF("Min time \t= %lf ms\n", min_time*1000);
	PRINTF("Max time \t= %lf ms\n", max_time*1000);
	PRINTF("Nb operation\t= %g Gflop\n", nb_operation);
	PRINTF("Performance\t= %g Gflop/s\n", nb_operation/(total/nbRUN));
	PRINTF("Error\t\t= %g\n", error/sum);
	//ADDCSV("AVG;%g;GFLOPS;%g;ERROR;%g\n",(total/nbRUN)*1000,nb_operation/(total/nbRUN),error/sum);

	/*************************************************************************************************************************/

	/* Liberation de la memoire **********************************************************************************************/
	FREE(signal_in);
	FREE(signal_out);
	FREE(reference);
	/*************************************************************************************************************************/
}
}

void readData(MKL_Complex8 *signal, size_t size_sig, MKL_Complex8 *reference, size_t size_ref) {
	MKL_Complex16 *tmp = (MKL_Complex16*) malloc(size_ref*sizeof(MKL_Complex16));
	if(tmp == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

    // a mettre a jour .....

	char *path = "../data/CovMatEstimation/CME_Complex_out_300.bin";

	FILE *fid = fopen(path, "rb");
	if(fid == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	size_t count = fread(tmp, sizeof(MKL_Complex16), size_ref, fid);
	if(count < (size_t) size_ref) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	for(size_t i = 0 ; i < (size_t) size_ref ; i++) {
		(*(reference+i)).real = (float) (*(tmp+i)).real;
		(*(reference+i)).imag = (float) (*(tmp+i)).imag;
	}

	free(tmp);
	fclose(fid);

	tmp = (MKL_Complex16*) malloc(size_sig*sizeof(MKL_Complex16));
	if(tmp == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	path = "../data/CovMatEstimation/CME_Complex_in_300.bin";

	fid = fopen(path, "rb");
	if(fid == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	count = fread(tmp , sizeof(MKL_Complex16), size_sig, fid);
	if(count < (size_t) size_sig) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}

	for(size_t i = 0 ; i < (size_t) size_sig ; i++) {
		(*(signal+i)).real = (float) (*(tmp+i)).real;
		(*(signal+i)).imag = (float) (*(tmp+i)).imag;
	}

	fclose(fid);
}
