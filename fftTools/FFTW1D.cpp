
#include "FFTW1D.H"
#include "fftw3.h"
#include <cstdlib>

using namespace std;

FFTW1D::FFTW1D(unsigned int a_N): FFT1D(a_N)
{
  m_in.resize(m_N);
  m_out.resize(m_N);
  fftw_complex* in; fftw_complex* out;
  in  = reinterpret_cast<fftw_complex*>(&(m_in[0]));
  out = reinterpret_cast<fftw_complex*>(&(m_out[0]));

  m_forward= fftw_plan_dft_1d(m_N, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
  m_inverse= fftw_plan_dft_1d(m_N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

FFTW1D::~FFTW1D()
{
  fftw_destroy_plan(m_forward);
  fftw_destroy_plan(m_inverse);
}

void FFTW1D::forwardFFTCC(vector<complex<double> > & a_fHat, 
                            const vector<complex<double> >& f) const 
{
  m_in = f;
  fftw_execute(m_forward);
  a_fHat = m_out;

}
  // inverse FFT: a_f[j] = \sum_{k=0}^{N-1} a_fHat[k] z^{j k}, z = e^{2 \pi \iota /m_N}
void FFTW1D::inverseFFTCC(vector<complex<double> > & a_f, 
                          const vector<complex<double> > & a_fHat) const
{
  m_in = a_fHat;
  fftw_execute(m_inverse);
  a_f = m_out;

}
