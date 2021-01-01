#include <stdio.h>

#include "benchmarks.h"

void display_benchmark_results(benchmark_results r) {
    char *used_precision = (r.precision == SP) ? "Simple Precsion" : "Double Precision";

    printf("*************\t%s %s\t*************\n", r.name, used_precision);
    printf("Params\t\t: %s\n", r.extra_infos);
	printf("Avg time\t= %lf ms\n", r.avg_time * 1000);
	printf("Min time \t= %lf ms\n", r.min_time*1000);
	printf("Max time \t= %lf ms\n", r.max_time*1000);
	printf("Nb operations\t= %i Gflop\n", r.gflop);
	printf("Performance\t= %g Gflop/s\n", r.gflops);
}