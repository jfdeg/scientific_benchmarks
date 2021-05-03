#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cufft.h>

#include "gene_bruit_rayleigh_scalaire.c"


int main (int argc, char** argv)
{

cudaSetDevice(0);

int NFFT,NPOINTS,NTOT;

NTOT=67108864;
//NTOT=16777216;


/* Chrono */
struct timespec tpdeb,tpfin;
clockid_t clock_id=CLOCK_REALTIME;
int status;
float dureeloc;
float dureetot = 0.0;
cufftResult_t res;

float time1,time2; 
cudaEvent_t start1,start2, stop1,stop2;


FILE *fichier1;
fichier1=fopen("../results/GPU_GFLOPS.dat","w");

/* Generateur de signal aleatoire Gaussien */
float param;
float *scalar_real,*scalar_imag;
float mflops,mflops2;
int i,j;
param=20.0;

//cudaSetDevice(1);

cufftComplex *data,*data_device,*data_device2;
cudaMalloc((void**)&data_device, sizeof(cufftComplex)*NTOT);
cudaMalloc((void**)&data_device2, sizeof(cufftComplex)*NTOT);

data=(cufftComplex *)calloc(NTOT,sizeof(cufftComplex));
scalar_real=(float *)calloc(NTOT,sizeof(float));
scalar_imag=(float *)calloc(NTOT,sizeof(float));


 // Debut boucle 
printf("Generation du bruit\n");
for(j = 0 ; j<NTOT ; j++)
{
gene_bruit_rayleigh_scalaire(param,scalar_real+i,scalar_imag+i);
}

for(i = 0 ; i < NTOT ; i++ ) 
	{
        data[i].x = scalar_real[i];
        data[i].y = scalar_imag[i];
        }
for (NPOINTS=2 ; NPOINTS<262144+1 ; NPOINTS=2*NPOINTS)
{
NFFT=NTOT/NPOINTS;
/* Declaration du plan et des donnes */
cufftHandle plan;


	/*  Creation du plan cufft */

	res=cufftPlan1d(&plan,NPOINTS, CUFFT_C2C, NFFT);
        printf("Allocation du plan cuFFT = %i\n",res);
	// Affichage du temps


/*
for (i = 0 ; i<NPOINTS ; i++)
{
fprintf(fichier1,"%20.15e\n",data[i].x);
fprintf(fichier2,"%20.15e\n",data[i].y);
}
*/
		cudaEventCreate(&start1);
		cudaEventCreate(&start2);
		cudaEventCreate(&stop1);
		cudaEventCreate(&stop2);
		cudaEventSynchronize( start1 );
		cudaEventRecord(start1, 0);


	cudaMemcpy(data_device, data, sizeof(cufftComplex)*NPOINTS*NFFT,cudaMemcpyHostToDevice);

        status=clock_gettime(clock_id, &tpdeb);
	cudaEventSynchronize( start2 );
	cudaEventRecord(start2, 0);
	/*----------- Debut calcul de la FFT----------*/

	/* Copie dans la memoire GPU */
	cudaThreadSynchronize();

	res = cufftExecC2C(plan,data_device,data_device2,CUFFT_FORWARD);
//	cufftExecC2C(plan,data_device,data_device,CUFFT_INVERSE);

	cudaThreadSynchronize();

	cudaEventSynchronize( stop2 );
	cudaEventRecord( stop2, 0);

	cudaEventSynchronize( stop2 );

	/*---------------------------------------------*/

	cudaMemcpy(data, data_device, sizeof(cufftComplex)*NPOINTS*NFFT,cudaMemcpyDeviceToHost);
	cudaThreadSynchronize();
        status=clock_gettime(clock_id, &tpfin);

	cudaEventSynchronize( stop1 );
	cudaEventRecord( stop1, 0);
	cudaEventElapsedTime( &time1, start1, stop1 );
	cudaEventElapsedTime( &time2, start2, stop2 );
	cudaEventDestroy( start1 );
	cudaEventDestroy( start2 );
	cudaEventDestroy( stop1 );
	cudaEventDestroy( stop2 );

        if (status<0)  fprintf(stderr,"Erreur clock_gettime (f:%s n:%d)\n",__FILE__,__LINE__);
        if (status<0)   printf("Erreur CLOCKGETTIME 2");
        dureeloc=(float)(tpfin.tv_sec-tpdeb.tv_sec)+(float)(tpfin.tv_nsec-tpdeb.tv_nsec)*1.e-9;
	dureetot=dureetot+dureeloc;

        printf("Execution de cuFFT = %i\n",res);
	// Affichage du temps

printf("Temps pour %d FFT de %d points complexes = %f ms \n",NFFT,NPOINTS,time2);
printf("Temps GPU avec transfert (timeur CPU) = %f ms \n",dureeloc*1000);
time2=time2/1000;
mflops=NFFT*(5*NPOINTS*(log10(NPOINTS)/log10(2)))/dureeloc;
mflops2=NFFT*(5*NPOINTS*(log10(NPOINTS)/log10(2)))/time2;
printf("Puissance sans transfert = %f Gflops \n",mflops2/1e9);
fprintf(fichier1,"%20.15f\n",mflops2/1e9);


cufftDestroy(plan);



}
free(data);
free(scalar_real);
free(scalar_imag);
cudaFree(data_device);
fclose(fichier1);

return 0;
}
