
#include <fstream>
#include <iterator>
#include <algorithm>

#include "gtest/gtest.h"

extern "C" {
#include "utils.h"
}

TEST(CovMatMKL, basic) {
    std::string dirpath = "tests/resources/data/";

    std::array files{
        "testcovmat_sp_in_10_10_real.bin",
        "testcovmat_sp_in_10_10_imag.bin",
        "testcovmat_sp_out_10_10_real.bin",
        "testcovmat_sp_out_10_10_imag.bin"
    };

    uint n_file = 0;
    uint const size = 10;
    uint const n_elem = 100;

    float data_ref[4][n_elem];

    for (auto file : files) {
        std::ifstream in_file(dirpath + file, std::ios::in | std::ios::binary);

        for (uint i = 0; i < n_elem; i++) {
            in_file.read((char*) data_ref[n_file], n_elem * sizeof(float));
        }

        in_file.close();
        n_file++;
    }

    ref_dataset_sp dataset = read_data_sp((char*) "tests/resources/data", (char*) "testcovmat", size, size);

    ASSERT_EQ(dataset.n, size);
    ASSERT_EQ(dataset.m, size);

    for (uint i = 0; i < n_elem; i++) {
        ASSERT_EQ(dataset.in_data_real[i], data_ref[0][i]);
        ASSERT_EQ(dataset.in_data_imag[i], data_ref[1][i]);
        ASSERT_EQ(dataset.out_data_real[i],data_ref[2][i]);
        ASSERT_EQ(dataset.out_data_imag[i],data_ref[3][i]);
    }
}
