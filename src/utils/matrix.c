#include <stdlib.h>

#include "matrix.h"

ref_dataset_sp new_ref_dataset_sp(size_t n, size_t m) {
    ref_dataset_sp dataset;
    
    dataset.n = n;
    dataset.m = m;

    dataset.in_data_real = malloc(n * m * sizeof(float));
    dataset.in_data_imag = malloc(n * m * sizeof(float));
    dataset.out_data_real = malloc(n * m * sizeof(float));
    dataset.out_data_imag = malloc(n * m * sizeof(float));

    return dataset;
}

void free_ref_dataset_sp(ref_dataset_sp dataset) {
    free(dataset.in_data_real);
    free(dataset.in_data_imag);
    free(dataset.out_data_real);
    free(dataset.out_data_imag);

    return;
}
