import numpy as np

import covmat_test_data

SIZES = [10]

if __name__ == "__main__":

    for size in SIZES:
        covmat_test_data.generate(size)