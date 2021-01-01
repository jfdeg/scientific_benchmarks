#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cublas_v2.h"
#include <cuda_runtime.h>
#include "gene_bruit_rayleigh_scalaire.c"


int main (void)
{


int M,N;
int kboucle;

cudaSetDevice(1);

FILE *fichier1;

fichier1=fopen("../results/results_cuda_gemm_fp32.dat","w");

/* Debut grand boucle */

for (kboucle=1 ; kboucle<31 ;kboucle++)
{

M=100*kboucle;
N=M;

printf(">>>>> Matrix size %dx%d  <<<<<<\n",M,N);

int iboucle,dim;
int param=20;
float *mat_real,*mat_imag;
float charge;

cublasStatus_t status;

// Chronometre

  struct timespec tpdeb,tpfin;
  clockid_t clock_id=CLOCK_REALTIME;
  int status2;

  struct timespec tpdeb2,tpfin2;

  float dureeloc,dureetot;
  dureetot=0.0;

// BLAS

cublasOperation_t transa,transb;


  /* CUBLAS */

float time1; 
cudaEvent_t start1, stop;

cuComplex cualpha,cubeta;
cuComplex *h_A,*h_B;
cuComplex *h_C;
cuComplex* d_A;
cuComplex* d_B;
cuComplex* d_C;

transa=CUBLAS_OP_N;
transb=CUBLAS_OP_N;


cualpha.x=1.0;
cualpha.y=0.0;
cubeta.x=0.0;
cubeta.y=0.0;

dim=M;

mat_real=(float*)calloc(M*N,sizeof(float));
mat_imag=(float*)calloc(M*N,sizeof(float));

h_A=(cuComplex*)calloc(M*N,sizeof(cuComplex));
h_B=(cuComplex*)calloc(M*N,sizeof(cuComplex));
h_C=(cuComplex*)calloc(M*N,sizeof(cuComplex));

    /* Initialize CUBLAS */

cublasHandle_t handle;
status=cublasCreate(&handle);


 cudaMalloc((void**)&d_A, M*N*sizeof(cuComplex));
 cudaMalloc((void**)&d_B, M*N*sizeof(cuComplex));
 cudaMalloc((void**)&d_C, M*N*sizeof(cuComplex));


for (iboucle=0 ; iboucle<M*N ; iboucle++)
{
gene_bruit_rayleigh_scalaire(param,mat_real+iboucle,mat_imag+iboucle);
}


for(iboucle=0 ; iboucle<N*N ; iboucle++)
{
h_A[iboucle].x=mat_real[iboucle];
h_A[iboucle].y=mat_imag[iboucle];
h_B[iboucle].x=mat_imag[iboucle];
h_B[iboucle].y=mat_real[iboucle];
}


  // Remise a zero et deuxieme chorno



  dureeloc=0.0;
  dureetot=0.0;

  cudaEventCreate(&start1);
  cudaEventCreate(&stop);
  status2=clock_gettime(clock_id, &tpdeb);

  // Copie de la matrice dans le GPU

  status = cublasSetVector(M*N, sizeof(cuComplex), h_A, 1, d_A, 1);

  status2=clock_gettime(clock_id, &tpdeb2);
  cudaThreadSynchronize();
  cudaEventRecord(start1, 0);

  cublasCgemm(handle,transa,transb,dim,dim,dim,&cualpha,d_A,dim,d_B,dim,&cubeta,d_C,dim);


	cudaThreadSynchronize();
	cudaEventRecord( stop, 0);

  status2=clock_gettime(clock_id, &tpfin2);

   status = cublasGetVector(M*N, sizeof(cuComplex), d_C, 1, h_C, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf (stderr, "!!!! device access error (read C)\n");
        return EXIT_FAILURE;
    }

  status2=clock_gettime(clock_id, &tpfin);
  dureeloc=(float)(tpfin.tv_sec-tpdeb.tv_sec)+(float)(tpfin.tv_nsec-tpdeb.tv_nsec)*1.e-9;
  dureetot=dureetot+dureeloc;

		cudaEventElapsedTime( &time1, start1, stop );
		cudaEventDestroy( start1 );
		cudaEventDestroy( stop );

  charge=(float)M;
  charge=8*charge*charge*charge;
  charge=charge/(time1/1000);
  printf("compute power = %f GFLOPS\n",charge/1e9);

  fprintf(fichier1,"%5.2f\n",charge/1e9);

cudaFree(d_A);
cudaFree(d_B);
cudaFree(d_C);

cublasDestroy ( handle ) ;

free(mat_real);
free(mat_imag);

free(h_A);
free(h_B);
free(h_C);

}

fclose(fichier1);
return 0;
}
