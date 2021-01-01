#include <math.h>
#define FORCE_OPENBLAS_COMPLEX_STRUCT
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef int MKL_INT;
typedef openblas_complex_double MKL_Complex16;

#define AVX512_ALIGN 64

#include "gene_bruit_rayleigh_scalaire.c"

int main(void) 
{

  int     val;
  int 	  M, N, K;
  int 	  lda, ldb, ldc;
  int     kboucle;
  float   charge;
  char    *p;
  
  FILE *fichier1;
  fichier1 = fopen("../results/results_openblas_gemm_fp64.dat", "w");

  int has_param_m = 0, has_param_n = 0, has_param_k = 0;
  char transa = CblasNoTrans, transb = CblasNoTrans;
  
  if ((p = getenv("MATMUL_M"))) {
    M = atoi(p);
    has_param_m=1;
  }
  if ((p = getenv("MATMUL_N"))) {
    N = atoi(p);
    has_param_n=1;
  } 
  if ((p = getenv("MATMUL_K"))) {
    K = atoi(p);
    has_param_k=1;
  }
  
  if ((p = getenv("TRANSA"))) {
    if(*p == 'T')
		transa = CblasTrans;
  }
  if ((p = getenv("TRANSB"))) {
    if(*p == 'T')
		transb = CblasTrans;
  }

  /* Begin huge loop */
  for (kboucle = 1; kboucle < 31; kboucle++) {

    val = 100 * kboucle;
    M = has_param_m?M:val;
    N = has_param_n?N:val;
	K = has_param_k?K:val;
	
	lda = (transa == CblasNoTrans)?M:K;
    ldb = (transb == CblasNoTrans)?K:N;
    ldc = M;
	
    printf(">>>>> Matrix size %dx%d (K=%d, TransA=%d TransB=%d) <<<<<<\n", M, N, K,\
		transa - CblasNoTrans, transb - CblasNoTrans);

    int    iboucle, jboucle;
    int    param = 20;
    float *mat_real, *mat_imag;

    // Timers
    struct timespec tpdeb, tpfin, tpcour;
    clockid_t clock_id = CLOCK_REALTIME;
    int status2;

    double dureeloc, dureetot;
    dureetot = 0.0;

    // BLAS
    MKL_Complex16 alpha, beta;
    MKL_Complex16 *A, *B, *C;

    alpha.real = 1.0;
    alpha.imag = 0.0;
    beta.real  = 0.0;
    beta.imag  = 0.0;

    mat_real = (float *)aligned_alloc(AVX512_ALIGN, val * val * sizeof(float));
    mat_imag = (float *)aligned_alloc(AVX512_ALIGN, val * val * sizeof(float));

    A = (MKL_Complex16 *)aligned_alloc(AVX512_ALIGN, M * K * sizeof(MKL_Complex16));
    B = (MKL_Complex16 *)aligned_alloc(AVX512_ALIGN, K * N * sizeof(MKL_Complex16));
    C = (MKL_Complex16 *)aligned_alloc(AVX512_ALIGN, M * N * sizeof(MKL_Complex16));

    for (iboucle = 0; iboucle < val * val; iboucle++) {
      gene_bruit_rayleigh_scalaire(param, mat_real + iboucle,
                                   mat_imag + iboucle);
    }

    // Force hermitian matrix
    for (iboucle = 0; iboucle < K; iboucle++) {
      for (jboucle = 0; jboucle < N; jboucle++) {
        B[(iboucle * N) + jboucle].real = (double)(*(mat_real + (jboucle * K) + iboucle));

        B[(iboucle * N) + jboucle].imag =
            -1.0 * (double)(*(mat_imag + (jboucle * K) + iboucle));
      }
    }

    for (iboucle = 0; iboucle < M * K; iboucle++) {
      A[iboucle].real = (double)mat_real[iboucle];
      A[iboucle].imag = (double)mat_imag[iboucle];
    }

    status2 = clock_gettime(clock_id, &tpdeb); // start timer;

    cblas_zgemm(CblasColMajor, transa, transb, \
        M, N, K, \
        &alpha, A, lda, \
        B, ldb, &beta, C, ldc);

    status2 = clock_gettime(clock_id, &tpfin); // stop timer

    dureeloc = (float)(tpfin.tv_sec - tpdeb.tv_sec) +
               (float)(tpfin.tv_nsec - tpdeb.tv_nsec) * 1.e-9;
    dureetot = dureetot + dureeloc;

    printf("execution time = %f ms\n", dureeloc * 1000);
    charge = (float)M;
    charge = 8.0f * charge * charge * charge;
    charge = charge / dureeloc;
    printf("%f GFLOPS\n", charge / 1e9);
    fprintf(fichier1, "%5.2f\n", charge / 1e9);

    free(mat_real);
    free(mat_imag);

    free(A);
    free(B);
    free(C);
  }
  fclose(fichier1);

  return 0;
}
