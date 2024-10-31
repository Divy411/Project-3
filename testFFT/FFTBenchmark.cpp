#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <memory>
#include <sys/time.h>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFT1DBRI.H"
#include "FFT1DW.H" 
#include "RectMDArray.H"
#include "DBox.H"
#include "FFTMD.H"



using namespace std;

double read() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if (!initialized) {
        gettimeofday(&start, NULL);
        initialized = true;
    }
    gettimeofday(&end, NULL);
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

// Copy of the `test` function from FFTMDTest.cpp
double test_new(const FFTMD& a_fftmd, double& a_seconds) {
    int N = a_fftmd.getN();
    int low[DIM], high[DIM];

    for (int dir = 0; dir < DIM; dir++) {
        low[dir] = 0;
        high[dir] = N - 1;
    }
    DBox b(low, high);

    RectMDArray<complex<double>> f(b);
    RectMDArray<complex<double>> fSave(b);
    RectMDArray<complex<double>> fHat(b);
    double h = 1.0 / N;

    for (Point pt = b.getLowCorner(); b.notDone(pt); b.increment(pt)) {
        double x = 1.0;
        for (int dir = 0; dir < DIM; dir++) {
            double y = (pt[dir] * h - 0.5);
            x += y * y * 32 * 32;
        }
        f[pt] = complex<double>(exp(-x), 0);
        fSave[pt] = f[pt];
    }

    double maxramp = 0.0;
    double maxiamp = 0.0;
    double maxamp = 0.0;
    Point ptmax = getZeros();
    double time = read();
    a_fftmd.forwardCC(f);
    a_seconds += read() - time;
    for (Point pt = b.getLowCorner(); b.notDone(pt); b.increment(pt)) {
        if (real(f[pt]) * real(f[pt]) + imag(f[pt]) * imag(f[pt]) > maxamp) {
            maxramp = real(f[pt]);
            maxiamp = imag(f[pt]);
            maxamp = real(f[pt]) * real(f[pt]) + imag(f[pt]) * imag(f[pt]);
            ptmax = pt;
        }
    }
    time = read();
    a_fftmd.inverseCC(f);
    a_seconds += read() - time;
    double maxerr = 0.0;
    int normalize = Power(N, DIM);

    for (Point pt = b.getLowCorner(); b.notDone(pt); b.increment(pt)) {
        double err = fabs(real(f[pt]) / normalize - real(fSave[pt]));
        if (err > maxerr) {
            maxerr = err;
        }
    }
    return maxerr;
}

// Rest of your benchmark code
int main() {
    // Define your FFTMD benchmark for M = [4, 5, 6, 7, 8, 9] with FFTW and BRI 1D solvers
    int Ms[] = {4, 5, 6, 7, 8, 9};
    for (const int M : Ms) {
        cout << "Running benchmark for M = " << M << " with FFTW" << endl;
        shared_ptr<FFT1D> p_fftFFTW = make_shared<FFT1DW>(M);
        FFTMD fftmdFFTW(p_fftFFTW);
        double timeFFTW = 0.0;
        double maxErrorFFTW = test_new(fftmdFFTW, timeFFTW);
        cout << "FFTW: time = " << timeFFTW << " seconds, error = " << maxErrorFFTW << endl;

        cout << "Running benchmark for M = " << M << " with BRI" << endl;
        shared_ptr<FFT1D> p_fftBRI = make_shared<FFT1DBRI>(M);
        FFTMD fftmdBRI(p_fftBRI);
        double timeBRI = 0.0;
        double maxErrorBRI = test_new(fftmdBRI, timeBRI);
        cout << "BRI: time = " << timeBRI << " seconds, error = " << maxErrorBRI << endl;
    }
    return 0;
}
