#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <omp.h>
#include "gene_bruit_rayleigh_scalaire.c"

#ifdef MKL_LIB
#include "mkl.h"
#include "mkl_dfti.h"
#endif


int main(int argc, char** argv) {

int NFFT,NPOINTS,NTOT;
NTOT=67108864;

printf("-*/*-*/- Programme MKL/FFTw simple precision -*/*-*/-\n");

FILE *fichier1;
fichier1=fopen("../results/MKL_32_GFLOPS.dat","w");

for (NPOINTS=2 ; NPOINTS<262144+1 ; NPOINTS=2*NPOINTS)
{
NFFT=NTOT/NPOINTS;


    struct timespec tpdeb, tpfin, tpcour;
    clockid_t clock_id = CLOCK_REALTIME;
    int status2;
    float dureeloc;
    float dureetot = 0.0;
    float mflops;

    float param;
    float *scalar_real, *scalar_imag;
    int i, j, k;
    param = 20.0;


    float f_num = 30000;
    float f_0 = f_num / 4.0;



#ifdef MKL_LIB
    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand [8];
    MKL_Complex8 * data[8];
    MKL_Complex8 * dataCour = NULL;


    memset(hand, 0, 8 * sizeof (DFTI_DESCRIPTOR_HANDLE));
    memset(data, 0, 8 * sizeof (MKL_Complex8));

#endif
    scalar_real = (float *) calloc(NPOINTS, sizeof (float));
    scalar_imag = (float *) calloc(NPOINTS, sizeof (float));


    //FILE *fichier1,*fichier2;

    //fichier1=fopen("/home/jfd/PROGRAMMES/FFT/scalar_real.dat","w");
    //fichier2=fopen("/home/jfd/PROGRAMMES/FFT/scalar_imag.dat","w");

#ifndef MKL_LIB

    /*----------- Declaration talbeaux FFT -----------------*/
    // CPU #1
    fftwf_complex *data, *fft_result, *data2, *fft_result2, *data3, *fft_result3, *data4, *fft_result4;
    fftwf_plan plan_forward, plan_forward2, plan_forward3, plan_forward4;

    data = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    data2 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result2 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    data3 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result3 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    data4 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result4 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);

    plan_forward = fftwf_plan_dft_1d(NPOINTS, data, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward2 = fftwf_plan_dft_1d(NPOINTS, data2, fft_result2, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward3 = fftwf_plan_dft_1d(NPOINTS, data3, fft_result3, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward4 = fftwf_plan_dft_1d(NPOINTS, data4, fft_result4, FFTW_FORWARD, FFTW_ESTIMATE);

    // CPU #2
    fftwf_complex *data5, *fft_result5, *data6, *fft_result6, *data7, *fft_result7, *data8, *fft_result8;
    fftwf_plan plan_forward5, plan_forward6, plan_forward7, plan_forward8;

    data5 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result5 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    data6 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result6 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    data7 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result7 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    data8 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);
    fft_result8 = (fftwf_complex*) fftwf_malloc(sizeof (fftwf_complex) * NPOINTS);

    plan_forward5 = fftwf_plan_dft_1d(NPOINTS, data5, fft_result5, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward6 = fftwf_plan_dft_1d(NPOINTS, data6, fft_result6, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward7 = fftwf_plan_dft_1d(NPOINTS, data7, fft_result7, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_forward8 = fftwf_plan_dft_1d(NPOINTS, data8, fft_result8, FFTW_FORWARD, FFTW_ESTIMATE);
    /*------------------------------------------------------*/
#else
    for (i = 0; i < 8; i++) {
        if ((data[i] = (MKL_Complex8 *) mkl_malloc(NPOINTS *
                (size_t)sizeof (MKL_Complex8), 64)) == NULL) {
            exit(0);
        }
        status = DftiCreateDescriptor(&hand[i],
                DFTI_SINGLE,
                DFTI_COMPLEX,
                1,
                (MKL_LONG) NPOINTS);
        if (0 != status) {
            printf(" ERROR, status = %li\n", status);
            exit(0);
        }

        status = DftiCommitDescriptor(hand[i]);
        if (0 != status) {
            printf(" %d ERROR, status = %li\n", __LINE__, status);
            exit(0);
        }
    }

#endif

    // Debut boucle 


    for (i = 0; i < NPOINTS; i++) {
        gene_bruit_rayleigh_scalaire(param, scalar_real + i, scalar_imag + i);
        // scalar_real[i] = cos((2 * M_PI * i * f_0) / (f_num));
        //scalar_imag[i] = 0;
    }
    for (j = 0; j < 8; j++) {
        dataCour = data[j];
        for (i = 0; i < NPOINTS; i++) {
            dataCour[i].real = scalar_real[i];
            dataCour[i].imag = scalar_imag[i];
        }
#if 0
        data2[i][0] = scalar_real[i];
        data2[i][1] = scalar_imag[i];
        data3[i][0] = scalar_real[i];
        data3[i][1] = scalar_imag[i];
        data4[i][0] = scalar_real[i];
        data4[i][1] = scalar_imag[i];
        data5[i][0] = scalar_real[i];
        data5[i][1] = scalar_imag[i];
        data6[i][0] = scalar_real[i];
        data6[i][1] = scalar_imag[i];
        data7[i][0] = scalar_real[i];
        data7[i][1] = scalar_imag[i];
        data8[i][0] = scalar_real[i];
        data8[i][1] = scalar_imag[i];
#endif

    }
    /*
    for (i = 0 ; i<NPOINTS ; i++)
    {
    fprintf(fichier1,"%20.15e\n",data[i][0]);
    fprintf(fichier2,"%20.15e\n",data[i][1]);
    }
     */
    // Fin remplissage des talbeaux 

    // Debut du chrono
    status2 = clock_gettime(clock_id, &tpdeb);
    if (status2 < 0) fprintf(stderr, "Erreur clock_gettime (f:%s n:%d)\n", __FILE__, __LINE__);
    if (status2 < 0) printf("Erreur CLOCKGETIME");


    /*----------- Calcul de la FFT----------*/
#pragma omp parallel num_threads(NBCORES) default(shared) private(j)
#ifndef MKL_LIB
    {
#pragma omp sections 
        {
            // CPU #1
#pragma omp section 
            {
                for (j = 0; j < NFFT / 8; j++) {
                    fftwf_execute(plan_forward);
                }
            }
#pragma omp section 
            {
                for (j = NFFT / 8; j < NFFT / 4; j++) {
                    fftwf_execute(plan_forward2);
                }
            }
#pragma omp section 
            {
                for (j = NFFT / 4; j < (3 * NFFT / 8); j++) {
                    fftwf_execute(plan_forward3);
                }
            }
#pragma omp section 
            {
                for (j = (3 * NFFT / 8); j < NFFT / 2; j++) {
                    fftwf_execute(plan_forward4);
                }
            }
            // CPU #1
#pragma omp section 
            {
                for (j = NFFT / 2; j < (5 * NFFT / 8); j++) {
                    fftwf_execute(plan_forward5);
                }
            }
#pragma omp section 
            {
                for (j = (5 * NFFT / 8); j < (3 * NFFT / 4); j++) {
                    fftwf_execute(plan_forward6);
                }
            }
#pragma omp section 
            {
                for (j = (6 * NFFT / 8); j < (7 * NFFT / 8); j++) {
                    fftwf_execute(plan_forward7);
                }
            }
#pragma omp section 
            {
                for (j = (7 * NFFT / 8); j < NFFT; j++) {
                    fftwf_execute(plan_forward8);
                }
            }
        }
    }
    /*--------------------------------------*/
#else
#if 0
    for (j = 0; j < NFFT; j++) {
        status = DftiComputeForward(hand, data);
        if (0 != status) {
            printf(" %d ERROR, status = %li\n", __LINE__, status);
            exit(0);
        }
    }
#endif
#pragma omp sections 
    {
        // CPU #1
#pragma omp section 
        {
            int k = 0;
            for (j = 0; j < NFFT / 8; j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
#pragma omp section 
        {
            int k = 1;
            for (j = NFFT / 8; j < NFFT / 4; j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
#pragma omp section 
        {
            int k = 2;
            for (j = NFFT / 4; j < (3 * NFFT / 8); j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
#pragma omp section 
        {
            int k = 3;
            for (j = (3 * NFFT / 8); j < NFFT / 2; j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
        // CPU #1
#pragma omp section 
        {
            int k = 4;
            for (j = NFFT / 2; j < (5 * NFFT / 8); j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
#pragma omp section 
        {
            int k = 5;
            for (j = (5 * NFFT / 8); j < (3 * NFFT / 4); j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
#pragma omp section 
        {
            int k = 6;
            for (j = (6 * NFFT / 8); j < (7 * NFFT / 8); j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
#pragma omp section 
        {
            int k = 7;
            for (j = (7 * NFFT / 8); j < NFFT; j++) {
                status = DftiComputeForward(hand[k], data[k]);
                if (0 != status) {
                    printf(" %d ERROR, status = %li\n", __LINE__, status);
                    exit(0);
                }
            }
        }
    }

    /*--------------------------------------*/
#endif

    status2 = clock_gettime(clock_id, &tpfin);
    if (status2 < 0) fprintf(stderr, "Erreur clock_gettime (f:%s n:%d)\n", __FILE__, __LINE__);
    if (status2 < 0) printf("Erreur CLOCKGETTIME 2");
    dureeloc = (float) (tpfin.tv_sec - tpdeb.tv_sec)+(float) (tpfin.tv_nsec - tpdeb.tv_nsec)*1.e-9;
    dureetot = dureetot + dureeloc;


// Affichage du temps
printf("Temps pour %d FFT de %d points complexes = %f ms \n",NFFT,NPOINTS,(dureeloc*1000));
printf("Temps CPU sans transfert (timeur GPU) = %f s \n",dureeloc);

mflops=NFFT*(5*NPOINTS*(log10(NPOINTS)/log10(2)))/dureeloc;
printf("PUISSANCE DE CALCUL = %f Gflops \n",mflops/1e9);
fprintf(fichier1,"%20.15f\n",mflops/1e9);

printf("\n\n");

#ifndef MKL_LIB
    fftwf_destroy_plan(plan_forward);
    fftwf_destroy_plan(plan_forward2);
    fftwf_destroy_plan(plan_forward3);
    fftwf_destroy_plan(plan_forward4); // CPU1
    fftwf_free(data);
    fftwf_free(data2);
    fftwf_free(data3);
    fftwf_free(data4);
    fftwf_free(fft_result);
    fftwf_free(fft_result2);
    fftwf_free(fft_result3);
    fftwf_free(fft_result4);

    fftwf_destroy_plan(plan_forward5);
    fftwf_destroy_plan(plan_forward6);
    fftwf_destroy_plan(plan_forward7);
    fftwf_destroy_plan(plan_forward8); // CPU2
    fftwf_free(data5);
    fftwf_free(data6);
    fftwf_free(data7);
    fftwf_free(data8);
    fftwf_free(fft_result5);
    fftwf_free(fft_result6);
    fftwf_free(fft_result7);
    fftwf_free(fft_result8);
#else
    for (j = 0; j < 8; j++) {
        if (hand[j])
            DftiFreeDescriptor(&hand[j]);
        if (data[j])
            mkl_free(data[j]);
    }
#endif
    free(scalar_real);
    free(scalar_imag);
}
fclose(fichier1);
    return 0;
}

