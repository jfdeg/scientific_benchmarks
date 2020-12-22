#pragma once

#include <stdlib.h>

typedef struct ref_dataset_sp {
    size_t n;
    size_t m;
    float *in_data_real;
    float *in_data_imag;
    float *out_data_real;
    float *out_data_imag;
} ref_dataset_sp;

ref_dataset_sp new_ref_dataset_sp(size_t n, size_t m);
void free_ref_dataset_sp(ref_dataset_sp dataset);

