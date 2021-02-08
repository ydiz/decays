#pragma once

#include "kaon/kaon.h"

namespace Grid {

// For each x and tK, return amplitude as a function of v

void typeII_D1a(const std::vector<int> &x, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &pl, const LatticePropagator &ps, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagatorSite Lxx; // L(x, x)
  peekSite(Lxx, pl, x);

  vector<LatticePropagator> Fu(4, env.grid), Gv(4, env.grid); // F_mu(u, x)
  for(int mu=0; mu<4; ++mu) Fu[mu] = adj(pl) * gmu5[mu] * wl[tK];

  for(int nu=0; nu<4; ++nu) Gv[nu] = adj(ws[tK]) * gmu5[nu] * ps; 
  LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
  for(int nu=0; nu<4; ++nu) Gv[nu] = Gv[nu] * exp_factor;           // G_nu(v, x) *= exp(M_k * (v_0 - t_K))

  vector<LatticePropagator> Cv = conv_with_E_typeII(Fu, env.M_K, max_uv_sep);

  LatticePropagator A(env.grid); A = Zero();
  for(int nu=0; nu<4; ++nu) A += Cv[nu] * Gv[nu];
  A = A - adj(A);  // Contribution of K0_bar

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int rho=0; rho<4; ++rho) {
    rst_Q1 += trace(gL[rho] * Lxx) * trace(gL[rho] * A);
    rst_Q2 += trace(gL[rho] * Lxx * gL[rho] * A);
  }
}

void typeII_D1b(const std::vector<int> &x, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &pl, const LatticePropagator &ps, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagatorSite Lxx; // L(x, x)
  peekSite(Lxx, pl, x);

  vector<LatticePropagator> Fv(4, env.grid), Gu(4, env.grid); // F_mu(u, x)

  for(int nu=0; nu<4; ++nu) Fv[nu] = adj(pl) * gmu5[nu] * wl[tK];
  LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
  for(int nu=0; nu<4; ++nu) Fv[nu] = Fv[nu] * exp_factor;           // F_nu(v, x) *= exp(M_k * (v_0 - t_K))

  for(int mu=0; mu<4; ++mu) Gu[mu] = adj(ws[tK]) * gmu5[mu] * ps; 

  vector<LatticePropagator> Cv = conv_with_E_typeII(Gu, env.M_K, max_uv_sep);

  LatticePropagator A(env.grid); A = Zero();
  for(int nu=0; nu<4; ++nu) A += Fv[nu] * Cv[nu];
  A = A - adj(A);  // Contribution of K0_bar

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int rho=0; rho<4; ++rho) {
    rst_Q1 += trace(gL[rho] * Lxx) * trace(gL[rho] * A);
    rst_Q2 += trace(gL[rho] * Lxx * gL[rho] * A);
  }

}



}
