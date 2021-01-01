# Resources for testing benchmark


All test data follow this naming convention:

```sh
BENCHMARK_PRECISON_IO_N_M_TYPE.bin
```

e.g. `covmat_sp_out_100_100_real.bin`

where:

- BENCHMARK: the name of the benchmark to test against
- PRECISION: float precision, either `sp` or `dp`
- IO: either `input` data for the test or the expected `output` of the benchmark method
- N: number of rows of the matrix
- M: number of columns of the matrix
- TYPE: either real or imag