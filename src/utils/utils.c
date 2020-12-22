#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "utils.h"
#include "matrix.h"


void read_datafile_sp(char *filepath, float *data, size_t data_size) {
	FILE *fid = fopen(filepath, "rb");

	if (fid == NULL)

		fprintf(
			stderr, "[Error] could not open data file %s at %s:%d : %s\n",
			filepath, __FILE__, __LINE__, strerror(errno)
		);

	size_t count = fread(data, sizeof(float), data_size, fid);

	if (count != data_size) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}
}

void read_datafile_dp(char *filepath, double *data, size_t data_size) {
	FILE *fid = fopen(filepath, "rb");

	if (fid == NULL)
		fprintf(
			stderr, "[Error] could not open data file %s at %s:%d : %s\n",
			filepath, __FILE__, __LINE__, strerror(errno)
		);

	size_t count = fread(data, sizeof(double), data_size, fid);

	if (count < data_size) {
		fprintf(stderr, "Error at %s:%d : %s\n", __FILE__, __LINE__, strerror(errno));
	}
}

ref_dataset_sp read_data_sp(char *dirpath, char *benchmark_name, int n, int m) {
	ref_dataset_sp dataset = new_ref_dataset_sp(n, m);

	char in_data_real_filepath_format[256], in_data_real_filepath[256];
	char in_data_imag_filepath_format[256], in_data_imag_filepath[256];
	char out_data_real_filepath_format[256], out_data_real_filepath[256];
	char out_data_imag_filepath_format[256], out_data_imag_filepath[256];

	strcpy(in_data_real_filepath_format, dirpath);
	strcpy(in_data_imag_filepath_format, dirpath);
	strcpy(out_data_real_filepath_format, dirpath);
	strcpy(out_data_imag_filepath_format, dirpath);

	strcat(in_data_real_filepath_format, "/%s_sp_in_%u_%u_real.bin");
	strcat(in_data_imag_filepath_format, "/%s_sp_in_%u_%u_imag.bin");
	strcat(out_data_real_filepath_format, "/%s_sp_out_%u_%u_real.bin");
	strcat(out_data_imag_filepath_format, "/%s_sp_out_%u_%u_imag.bin");


	sprintf(in_data_real_filepath, in_data_real_filepath_format, benchmark_name, n, m);
	sprintf(in_data_imag_filepath, in_data_imag_filepath_format, benchmark_name, n, m);
	sprintf(out_data_real_filepath, out_data_real_filepath_format, benchmark_name, n, m);
	sprintf(out_data_imag_filepath, out_data_imag_filepath_format, benchmark_name, n, m);

	read_datafile_sp(in_data_real_filepath, dataset.in_data_real, n * m);
	read_datafile_sp(in_data_imag_filepath, dataset.in_data_imag, n * m);
	read_datafile_sp(out_data_real_filepath, dataset.out_data_real, n * m);
	read_datafile_sp(out_data_imag_filepath, dataset.out_data_imag, n * m);

	return dataset;
}