#pragma once

#include "kaon/kaon.h"

namespace Grid {

std::vector<LatticePropagator> conv_with_E_typeII(const std::vector<LatticePropagator> &F, double M_K, int max_uv_sep) {
  using namespace std;
  GridBase *grid = F[0].Grid();

  FFT theFFT((GridCartesian *)grid);

  // static vector<vector<LatticeComplex>> Euv_fft_conj(4, vector<LatticeComplex>(4, grid));  // // Euv_fft = (Fourier{Euv})^*, where Euv is EM factor 
  static vector<vector<LatticeComplex>> Euv_fft_conj;  // Euv_fft = (Fourier{Euv})^*, where Euv is EM factor 
  static bool initialized = false;
  if(!initialized) {
    Euv_fft_conj = calc_Euv_fft_conj(grid, M_K, max_uv_sep);
    initialized = true;
  }

  vector<LatticePropagator> F_fft(3, grid); // F_fft[3] is not used; so I do not calculate it
  for(int mu=0; mu<3; ++mu) theFFT.FFT_all_dim(F_fft[mu], F[mu], FFT::forward);

  vector<LatticePropagator> rst(4, grid);
  for(int mu=0; mu<4; ++mu) rst[mu] = Zero();

  for(int nu=0; nu<3; ++nu) {   // Here, nu != 3, because E_{3, mu} = 0
    for(int mu=0; mu<3; ++mu) { // Here, mu != 3, because E_{nu, 3} = 0
      if(mu==nu) continue;  // if mu==nu, E_{mu,mu} = 0
      rst[nu] += Euv_fft_conj[mu][nu] * F_fft[mu];
    }
  }

  for(int nu=0; nu<3; ++nu) theFFT.FFT_all_dim(rst[nu], rst[nu], FFT::backward); // rst[3] is always 0

  return rst;
}




}
