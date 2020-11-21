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

#ifdef OPEN
#include <cblas.h>
typedef openblas_complex_float MKL_Complex8;
typedef openblas_complex_double MKL_Complex16;
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

extern int mklEstimation_STAP_sp(MKL_Complex8 *signal_out, MKL_Complex8 *signal_in, int nbCD, int nbRV, int nbR, int nbV, int nbREF, int nbREC, int nbFFT);
extern void estimation_alloc_once_sp(int nbRV, int nbR, int nbV, int nbREF);
extern void estimation_free_once_sp(void);

static void readData(MKL_Complex8 *signal, size_t size_sig, MKL_Complex8 *reference, size_t size_ref);

void AutoBenchEstimation_MKL_sp() {

    for(int N = 100 ; i <= 3000 ; i = i + 100) {

	/** Data init **/
	size_t size_sig = (size_t) (N*N*nbRUN);
	size_t size_ref = (size_t) (N*N*nbRUN);

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

	/* Mise en place du bench MKL ********************************************************************************************/
	float total = 0.0f;
	float min_time = 65536.0f;
	float max_time = 0.0f;
	float elapsedTime  = 0.0f;
	float nb_operation = (float) (1e-9*N*N*8.0f);

	struct timespec tpdeb;
	struct timespec tpfin;
	clockid_t clock_id = CLOCK_REALTIME;

	FILE *fid = fopen("MatrixMul_MKL_sp.txt", "w+");
	if(fid == NULL) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
		return;
	}

	// Warmup
	memset(signal_out, 0, size_ref*sizeof(MKL_Complex8));
	mklEstimation_STAP_sp(signal_out, signal_in, N);


	for(int i = 0 ; i < nb_fois ; i++) {
		memset(signal_out, 0, size_ref*sizeof(MKL_Complex8));

		clock_gettime(clock_id, &tpdeb);

		mklMatrixMul_sp(signal_out, signal_in, N);

		clock_gettime(clock_id, &tpfin);

		elapsedTime = (float) (tpfin.tv_sec - tpdeb.tv_sec) + (float) ((float) (tpfin.tv_nsec - tpdeb.tv_nsec) * 1e-9);
		total += elapsedTime;
		max_time = fmaxf(max_time, elapsedTime);
		min_time = fminf(min_time, elapsedTime);
		fprintf(fid, "%g\n", elapsedTime*1000);
	}

	fclose(fid);



	/*************************************************************************************************************************/

	/* Comparaison des resultats MKL *****************************************************************************************/
	MKL_Complex8 diff;
	float error = 0.0f;
	float sum   = 0.0f;

	for(int i = 0 ; i < nbCD ; i++) {
		for(int j = 0 ; j < nbRV ; j++) {
			for(int k = 0 ; k < nbRV ; k++) {
				if(j >= k) {
					MKL_Complex8 *ref = (MKL_Complex8*) reference+k+(j+i*nbRV)*nbRV;
					MKL_Complex8 *sig = (MKL_Complex8*) signal_out+k+(j+i*nbRV)*nbRV;
					diff.real  = ref->real - sig->real;
					diff.imag  = ref->imag - sig->imag;
					error += sqrtf(diff.real * diff.real + diff.imag * diff.imag);
					sum   += sqrtf(ref->real * ref->real + ref->imag * ref->imag);
				}
			}
		}
	}

	PRINTF("*************\tTotal MKL Simple Precision\t*************\n");
	PRINTF("Average elapsed time\t= %lf ms\n", (total/nb_fois)*1000);
	PRINTF("Min time \t= %lf ms\n", min_time*1000);
	PRINTF("Max time \t= %lf ms\n", max_time*1000);
	PRINTF("Nb operation\t= %g Gflop\n", nb_operation);
	PRINTF("Performance\t= %g Gflop/s\n", nb_operation/(total/nb_fois));
	PRINTF("Error\t\t= %g\n", error/sum);
	ADDCSV("AVG;%g;GFLOPS;%g;ERROR;%g\n",(total/nb_fois)*1000,nb_operation/(total/nb_fois),error/sum);

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

	char *home = getenv("STAP_DATA_PATH");
	char path[150];
	char *type_fo = getenv("TYPE_FO");
	char path2[20];

	strcpy(path, home);
	strcpy(path2, type_fo);
	strcat(path, path2);
	strcat(path, "/estimation.bin");

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

	strcpy(path, home);
	strcpy(path2, type_fo);
	strcat(path, path2);
	strcat(path, "/reformation.bin");

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
