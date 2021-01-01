import numpy as np

np.random.seed(0)

def generate(size):
    (m_real, m_imag) = np.random.randn(size, size), np.random.randn(size, size)

    m_complex = m_real + 1j * m_imag

    m_covmat = np.cov(m_complex)

    output_file_name = "tests/resources/data/testcovmat_sp_out_{}_{}".format(size, size)
    input_file_name = "tests/resources/data/testcovmat_sp_in_{}_{}".format(size, size)

    m_real.astype("float32").tofile(input_file_name + "_real.bin")
    m_imag.astype("float32").tofile(input_file_name + "_imag.bin")

    np.real(m_covmat).astype("float32").tofile(output_file_name + "_real.bin")
    np.imag(m_covmat).astype("float32").tofile(output_file_name + "_imag.bin")
