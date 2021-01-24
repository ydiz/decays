#pragma once

namespace Grid {


// Calculate Euv_fft = (Fourier{Euv})^*, where Euv is EM factor on half lattice
std::vector<std::vector<LatticeComplex>> calc_Euv_fft_conj(GridBase *grid, double M_K, int max_uv_sep) {
    std::vector<std::vector<LatticeComplex>> Euv_fft_conj(4, std::vector<LatticeComplex>(4, grid));

    LatticeKGG Euv(grid);
    EM_factor_half_lattice(Euv, {0,0,0,0}, M_K, max_uv_sep);

    FFT theFFT((GridCartesian *)grid);
    for(int mu=0; mu<4; ++mu) { 
      for(int nu=0; nu<4; ++nu) { 
        if(mu==nu || mu==3 || nu == 3) Euv_fft_conj[mu][nu] = Zero(); // Under these conditions, E_{mu,nu} = 0
        else {
          LatticeComplex tmp = peekLorentz(Euv, mu, nu);
          theFFT.FFT_all_dim(Euv_fft_conj[mu][nu], tmp, FFT::forward);
          Euv_fft_conj[mu][nu] = conjugate(Euv_fft_conj[mu][nu]);
        }
      }
    }
    return Euv_fft_conj;
}





// return exp(M_K * (v_0 - tK)); if v_0 - tK > T/2, set to 0
LatticeComplex exp_v0_tK(GridCartesian *grid, int t_K, double M_K) {

  LatticeComplex lat(grid);
  int T = lat.Grid()->_fdimensions[3];

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

		LatticeComplexSite m;

    int left_dist = left_distance(gcoor[3], t_K, T); // v_0 - tK // the distance from v0 to tK if you can only move to the left
    if(left_dist <= T/2) {
      double val = std::exp(M_K * left_dist); 
      m()()() = Complex(val, 0.);
    }
    else m = Zero();

		pokeLocalSite(m, lat, lcoor);
  });

  return lat;
}








}
