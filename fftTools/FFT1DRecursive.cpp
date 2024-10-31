#include "FFT1DRecursive.H"
#include <vector>            // std::vector
#include <complex>           // std::complex
#include <iostream>          // std::cout, std::endl
#include <cmath>             // std::polar, M_PI
#include <cstdlib>           // std::abort

FFT1DRecursive::FFT1DRecursive(const unsigned int& a_M):FFT1D(a_M)
{
 // do you need to build any data structures to help implement this class?
} 


FFT1DRecursive::~FFT1DRecursive() 
{
  // be sure to clean up anything you created in the constructor.
  // most STL containers will clean themselves up for you.
}

void FFT1DRecursive::forwardFFTCC(std::vector<std::complex<double>>& a_fHat,
                                  const std::vector<std::complex<double>>& f) const {
    // Ensure input size matches expected size
    unsigned int N = f.size();
    a_fHat.resize(N);

    if (N <= 1) {
        // Base case: just copy the input when size is 1
        a_fHat = f;
        return;
    }

    // Divide: Separate even and odd elements
    std::vector<std::complex<double>> even(N / 2), odd(N / 2);
    for (unsigned int i = 0; i < N / 2; ++i) {
        even[i] = f[2 * i];
        odd[i] = f[2 * i + 1];
    }

    // Recursively solve for even and odd parts
    std::vector<std::complex<double>> evenFFT(N / 2), oddFFT(N / 2);
    forwardFFTCC(evenFFT, even);
    forwardFFTCC(oddFFT, odd);

    // Combine the results
    for (unsigned int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / N) * oddFFT[k];
        a_fHat[k] = evenFFT[k] + t;
        a_fHat[k + N / 2] = evenFFT[k] - t;
    }
}

void FFT1DRecursive::inverseFFTCC(std::vector<std::complex<double>>& a_f,
                                  const std::vector<std::complex<double>>& a_fHat) const {
    unsigned int N = a_fHat.size();
    a_f.resize(N);

    if (N <= 1) {
        a_f = a_fHat;
        return;
    }

    // Divide: Separate even and odd elements
    std::vector<std::complex<double>> even(N / 2), odd(N / 2);
    for (unsigned int i = 0; i < N / 2; ++i) {
        even[i] = a_fHat[2 * i];
        odd[i] = a_fHat[2 * i + 1];
    }

    // Recursively solve for even and odd parts
    std::vector<std::complex<double>> evenIFFT(N / 2), oddIFFT(N / 2);
    inverseFFTCC(evenIFFT, even);
    inverseFFTCC(oddIFFT, odd);

    // Combine the results
    for (unsigned int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar(1.0, 2 * M_PI * k / N) * oddIFFT[k];
        a_f[k] = (evenIFFT[k] + t) / 2.0;
        a_f[k + N / 2] = (evenIFFT[k] - t) / 2.0;
    }
}

