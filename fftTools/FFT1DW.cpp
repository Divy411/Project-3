#include <cstdlib>
#include "FFT1DW.H"

FFT1DW::FFT1DW(unsigned int a_N) : FFT1D(a_N) {
    // Initialize input and output arrays
    m_in.resize(m_N);
    m_out.resize(m_N);

    // Create FFTW plans for forward and inverse FFT
    fftw_complex* in = reinterpret_cast<fftw_complex*>(m_in.data());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(m_out.data());

    m_forward = fftw_plan_dft_1d(m_N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    m_inverse = fftw_plan_dft_1d(m_N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

FFT1DW::~FFT1DW() {
    // Clean up FFTW plans
    fftw_destroy_plan(m_forward);
    fftw_destroy_plan(m_inverse);
}

void FFT1DW::forwardFFTCC(std::vector<std::complex<double>>& a_fHat,
                          const std::vector<std::complex<double>>& f) const {
    // Copy input data into FFTW input array
    m_in = f;

    // Execute forward FFT
    fftw_execute(m_forward);

    // Copy result into output vector
    a_fHat = m_out;
}

void FFT1DW::inverseFFTCC(std::vector<std::complex<double>>& a_f,
                          const std::vector<std::complex<double>>& a_fHat) const {
    // Copy input data into FFTW input array
    m_in = a_fHat;

    // Execute inverse FFT
    fftw_execute(m_inverse);

    // Copy result into output vector
    a_f = m_out;
}
