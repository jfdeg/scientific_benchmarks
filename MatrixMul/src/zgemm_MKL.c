#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mkl_blas.h>

#include "gene_bruit_rayleigh_scalaire.c"



int main (void)
{

int M,N;
int kboucle;
float charge;


FILE *fichier1;
fichier1=fopen("../results/results_mkl_gemm_fp64.dat","w");

/* Debut grand boucle */

for (kboucle=1 ; kboucle<31 ;kboucle++)
{

M=100*kboucle;
N=M;

printf(">>>>> Matrix size %dx%d  <<<<<<\n",M,N);

int iboucle,jboucle;
int param=20;
float *mat_real,*mat_imag;

// Chronometre

  struct timespec tpdeb,tpfin,tpcour;
  clockid_t clock_id=CLOCK_REALTIME;
  int status2;

  float dureeloc,dureetot;
  dureetot=0.0;

// BLAS

char transa,transb;
MKL_INT dim;
MKL_INT lda;
MKL_Complex16 *alpha,*beta;
MKL_Complex16 *A,*B,*C;

transa='N';
transb='N';

alpha=(MKL_Complex16*)calloc(1,sizeof(MKL_Complex16));
beta=(MKL_Complex16*)calloc(1,sizeof(MKL_Complex16));

alpha[0].real=1.0;
alpha[0].imag=0.0;
beta[0].real=0.0;
beta[0].imag=0.0;

dim=M;

mat_real=(float*)calloc(M*N,sizeof(float));
mat_imag=(float*)calloc(M*N,sizeof(float));

A=(MKL_Complex16*)calloc(M*N,sizeof(MKL_Complex16));
B=(MKL_Complex16*)calloc(M*N,sizeof(MKL_Complex16));
C=(MKL_Complex16*)calloc(M*N,sizeof(MKL_Complex16));


for (iboucle=0 ; iboucle<M*N ; iboucle++)
{
gene_bruit_rayleigh_scalaire(param,mat_real+iboucle,mat_imag+iboucle);
}

// Force hermitian matrix

for (iboucle=0; iboucle<N; iboucle++)
  {
  for (jboucle=0; jboucle<N; jboucle++)
    {
    B[(iboucle*N)+jboucle].real=
    *(mat_real+(jboucle*N)+iboucle);

    B[(iboucle*N)+jboucle].imag=
    -1.*(*(mat_imag+(jboucle*N)+iboucle));
    }
  }

for(iboucle=0 ; iboucle<N*N ; iboucle++)
{
A[iboucle].real=(double)mat_real[iboucle];
A[iboucle].imag=(double)mat_imag[iboucle];

}


status2=clock_gettime(clock_id, &tpdeb); // start timer

zgemm(&transa,&transb,&dim,&dim,&dim,alpha,A,&dim,B,&dim,beta,C,&dim);

status2=clock_gettime(clock_id, &tpfin); // stop timer

dureeloc=(float)(tpfin.tv_sec-tpdeb.tv_sec)+(float)(tpfin.tv_nsec-tpdeb.tv_nsec)*1.e-9;
dureetot=dureetot+dureeloc;

  printf("execution time = %f ms\n",dureeloc*1000);
  charge=(float)M;
  charge=8*charge*charge*charge;
  charge=charge/dureeloc;
  printf("%f GFLOPS\n",charge/1e9);
  fprintf(fichier1,"%5.2f\n",charge/1e9);


free(mat_real);
free(mat_imag);

free(A);
free(B);
free(C);

}
fclose(fichier1);

return 0;
}
