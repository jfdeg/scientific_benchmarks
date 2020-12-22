#include <fstream>
#include <iterator>
#include <algorithm>

#include "gtest/gtest.h"

extern "C" {
#include "utils.h"
}

TEST(Utils, read_data_sp) {
    std::string dirpath = "tests/resources/data/";

    std::array files{
        "testread_sp_in_10_10_real.bin",
        "testread_sp_in_10_10_imag.bin",
        "testread_sp_out_10_10_real.bin",
        "testread_sp_out_10_10_imag.bin"
    };

    uint n_file = 0;
    uint const size = 10;
    uint const n_elem = 100;

    float f;

    for (auto file : files) {
        std::ofstream output(dirpath + file,  std::ios::out | std::ios::binary);

        for (uint i = 0; i < n_elem; i++) {
            f = (float) i + (n_elem * n_file);
            output.write((char*) &f, sizeof(float));
        }

        output.flush();
        output.close();
        n_file++;
    }

    ref_dataset_sp dataset = read_data_sp((char*) "tests/resources/data", (char*) "testread", size, size);

    ASSERT_EQ(dataset.n, size);
    ASSERT_EQ(dataset.m, size);

    for (uint i = 0; i < n_elem; i++) {
        ASSERT_EQ(dataset.in_data_real[i], (float) i);
        ASSERT_EQ(dataset.in_data_imag[i], (float) i + 100);
        ASSERT_EQ(dataset.out_data_real[i], (float) i + 200);
        ASSERT_EQ(dataset.out_data_imag[i], (float) i + 300);
    }
}
