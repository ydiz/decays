#pragma once

#include "kaon/kaon.h"

namespace Grid {

std::vector<std::vector<LatticeColourMatrix>> conv_with_E_typeI(const std::vector<std::vector<LatticeColourMatrix>> &F, double M_K, int max_uv_sep) {
  using namespace std;
  GridBase *grid = F[0][0].Grid();

  FFT theFFT((GridCartesian *)grid);

  static vector<vector<LatticeComplex>> Euv_fft_conj;  // Euv_fft = (Fourier{Euv})^*, where Euv is EM factor 
  static bool initialized = false;
  if(!initialized) {
    Euv_fft_conj = calc_Euv_fft_conj(grid, M_K, max_uv_sep);
    initialized = true;
  }

  vector<vector<LatticeColourMatrix>> F_fft(3, vector<LatticeColourMatrix>(4, grid)); // F_fft[3, rho] is not used; so I do not calculate it
  for(int mu=0; mu<3; ++mu) 
    for(int rho=0; rho<4; ++rho)
      theFFT.FFT_all_dim(F_fft[mu][rho], F[mu][rho], FFT::forward);

  vector<vector<LatticeColourMatrix>> rst(4, vector<LatticeColourMatrix>(4, grid)); 
  for(int nu=0; nu<4; ++nu) {   
    for(int rho=0; rho<4; ++rho) {

      rst[nu][rho] = Zero();

      if(nu==3) continue;   // rst[3][rho] is always 0

      for(int mu=0; mu<3; ++mu) {   // Here, mu != 3, because E_{nu, 3} = 0
        if(mu==nu) continue;  // if mu==nu, E_{mu,mu} = 0
        rst[nu][rho] += Euv_fft_conj[mu][nu] * F_fft[mu][rho];
      }
    }
  }

  for(int nu=0; nu<3; ++nu)  // rst[3][rho] is always 0
    for(int rho=0; rho<4; ++rho)
      theFFT.FFT_all_dim(rst[nu][rho], rst[nu][rho], FFT::backward);

  return rst;
}




}
