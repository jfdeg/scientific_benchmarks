#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <omp.h>
#include <mkl.h>


void mklEstimation_Complex_sp(MKL_Complex8 *signal_out, MKL_Complex8 *signal_in, int N, int M, int nbMat, char transN, char transC, MKL_Complex8 norm, MKL_Complex8 null) {
	
    for(int cpt = 0 ; cpt < nbMat ; cpt++){
        cgemm( &transN, &transC, &N, &N, &M, &norm, signal_in, &N, signal_in, &N, &null, signal_out, &N );
    }
}
