#pragma once

#include "matrix.h"

typedef enum verbosity {
    DEBUG,
    INFO,
    NONE
} verbosity;

ref_dataset_sp read_data_sp(char *dirpath, char *benchmark_name, int n, int m);
