#pragma once

#include "kaon/kaon.h"

namespace Grid {

void typeIV_D1a(const std::vector<int> &v, int tK, 
                LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, LatticeComplex &sBar_d_T2D1a,
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &pl, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticeKGG Euv(env.grid);
  EM_factor_half_lattice(Euv, v, env.M_K, max_uv_sep); // EM factor centered at v

  // g = \sum_u gnu g5 L(u, v)^dagger gmu g5 L(u, tK)
  LatticePropagator g_u(env.grid);  g_u = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      g_u += gmu5[nu] * adj(pl) * gmu5[mu] * wl[tK] * Euv_munu;
    }
  }
  LatticePropagatorSite g = sum(g_u);

  // calculate contraction
  LatticePropagator f_Q1_K(env.grid), f_Q1_Kbar(env.grid), f_Q2_K(env.grid), f_Q2_Kbar(env.grid), f_sBar_d(env.grid);
  f_Q1_K = Zero(); f_Q1_Kbar = Zero(); f_Q2_K = Zero(); f_Q2_Kbar = Zero();
  for(int rho=0; rho<4; ++rho) {
    f_Q1_K += trace(gL[rho] * Lxx) * adj(ws[tK]) * gL[rho] * pl;
    f_Q1_Kbar += trace(gL[rho] * Lxx) * adj(pl) * gL[rho] * ws[tK];

    f_Q2_K += adj(ws[tK]) * gL[rho] * Lxx *  gL[rho] * pl;
    f_Q2_Kbar += adj(pl) * gL[rho] * Lxx * gL[rho] * ws[tK];
  }
  f_sBar_d = adj(ws[tK]) * g5 * pl;

  rst_Q1 = trace(f_Q1_K * g + f_Q1_Kbar * adj(g));
  rst_Q2 = trace(f_Q2_K * g + f_Q2_Kbar * adj(g));
  sBar_d_T2D1a = trace(f_sBar_d * g); 

  int tsep = (v[3] - tK + T) % T;
  rst_Q1 *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
  rst_Q2 *= std::exp(env.M_K * tsep);
  sBar_d_T2D1a *= std::exp(env.M_K * tsep);
}



void typeIV_D2a(const std::vector<int> &v, int tK, 
                LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, LatticeComplex &sBar_d_T2D2a,
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &ps, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticeKGG Euv(env.grid);
  EM_factor_half_lattice(Euv, v, env.M_K, max_uv_sep); // EM factor centered at v

  // g = \sum_u gnu g5 L(u, v)^dagger gmu g5 L(u, tK)
  LatticePropagator g_u(env.grid);  g_u = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      g_u += adj(ws[tK]) * gmu5[mu] * ps * gmu5[nu] * Euv_munu;
    }
  }
  LatticePropagatorSite g = sum(g_u);

  // calculate contraction
  LatticePropagator f_Q1_K(env.grid), f_Q1_Kbar(env.grid), f_Q2_K(env.grid), f_Q2_Kbar(env.grid), f_sBar_d(env.grid);
  f_Q1_K = Zero(); f_Q1_Kbar = Zero(); f_Q2_K = Zero(); f_Q2_Kbar = Zero();
  for(int rho=0; rho<4; ++rho) {
    f_Q1_K += trace(gL[rho] * Lxx) * adj(ps) * gL[rho] * wl[tK];
    f_Q1_Kbar += trace(gL[rho] * Lxx) * adj(wl[tK]) * gL[rho] * ps;

    f_Q2_K += adj(ps) * gL[rho] * Lxx *  gL[rho] * wl[tK];
    f_Q2_Kbar += adj(wl[tK]) * gL[rho] * Lxx * gL[rho] * ps;
  }
  f_sBar_d = adj(ps) * g5 * wl[tK];

  rst_Q1 = trace(f_Q1_K * g + f_Q1_Kbar * adj(g));
  rst_Q2 = trace(f_Q2_K * g + f_Q2_Kbar * adj(g));
  sBar_d_T2D2a = trace(f_sBar_d * g); 

  int tsep = (v[3] - tK + T) % T;
  rst_Q1 *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
  rst_Q2 *= std::exp(env.M_K * tsep);
  sBar_d_T2D2a *= std::exp(env.M_K * tsep);
}





}

