#pragma once

#include "kaon/kaon.h"

namespace Grid {

void typeV(const std::vector<int> &v, int tK, LatticeComplex &rst_D1Q1, LatticeComplex &rst_D1Q2, 
           LatticeComplex &rst_D2Q1, LatticeComplex &rst_D2Q2, LatticeComplex &sBar_d_D1, LatticeComplex &sBar_d_D2,
           Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
           const LatticePropagator &pl, const LatticePropagator &ps, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticeKGG Euv(env.grid);
  EM_factor_half_lattice(Euv, v, env.M_K, max_uv_sep); // EM factor centered at v

  LatticeComplex fu_D1(env.grid), fu_D2(env.grid);
  fu_D1 = Zero(); fu_D2 = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      fu_D1 += trace(gmu5[mu] * pl * gmu5[nu] * adj(pl)) * Euv_munu;
      fu_D2 += trace(gmu5[mu] * ps * gmu5[nu] * adj(ps)) * Euv_munu;
    }
  }
  LatticeComplexSite f_D1 = sum(fu_D1);
  LatticeComplexSite f_D2 = sum(fu_D2);

  LatticePropagator g_x = wl[tK] * adj(ws[tK]); 

  // sBar_d_D1 = f_D1 * trace(g5 * g_x);
  // sBar_d_D2 = f_D2 * trace(g5 * g_x);
  sBar_d_D1 = f_D1 * 2. * real(trace(g5 * g_x)); // 2. * real(xxx) is for adding contribution of K0bar
  sBar_d_D2 = f_D2 * 2. * real(trace(g5 * g_x));

  g_x = g_x + adj(g_x);  // add contribution of K0bar

  rst_D1Q1 = Zero(); rst_D1Q2 = Zero(); rst_D2Q1 = Zero(); rst_D2Q2 = Zero();
  for(int rho=0; rho<4; ++rho) {
    rst_D1Q1 += f_D1 * trace(gL[rho] * g_x)  * trace(gL[rho] * Lxx);
    rst_D1Q2 += f_D1 * trace(gL[rho] * g_x * gL[rho] * Lxx);
    rst_D2Q1 += f_D2 * trace(gL[rho] * g_x)  * trace(gL[rho] * Lxx);
    rst_D2Q2 += f_D2 * trace(gL[rho] * g_x * gL[rho] * Lxx);
  }

  int tsep = (v[3] - tK + T) % T;
  sBar_d_D1 *= std::exp(env.M_K * tsep);
  sBar_d_D2 *= std::exp(env.M_K * tsep);
  rst_D1Q1 *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
  rst_D1Q2 *= std::exp(env.M_K * tsep);
  rst_D2Q1 *= std::exp(env.M_K * tsep); 
  rst_D2Q2 *= std::exp(env.M_K * tsep);
}

}
