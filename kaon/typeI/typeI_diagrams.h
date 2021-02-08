#pragma once

#include "kaon/kaon.h"

namespace Grid {

// For each x and tK, return amplitude as a function of v
// Possible optimization for type I: for diagram b, Fourier transformation does not depend on tK and can be done outside of the tsep loop

void typeI_D1a(const std::vector<int> &x, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
               Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
               const LatticePropagator &pl, const LatticePropagator &ps, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagatorSite ws_x; // H(x, tK)
  peekSite(ws_x, ws[tK], x);

  vector<vector<LatticeColourMatrix>> Fu(4, vector<LatticeColourMatrix>(4, env.grid)); // F_{mu,rho}(u, x)
  for(int mu=0; mu<4; ++mu) {
    LatticePropagator A = adj(pl) * gmu5[mu] * wl[tK] * adj(ws_x);
    A = A - adj(A);              //  add contribution from K0_bar
    for(int rho=0; rho<4; ++rho) Fu[mu][rho] = traceS(gL[rho] * A); 
  }

  vector<vector<LatticeColourMatrix>> Gv(4, vector<LatticeColourMatrix>(4, env.grid)); // G_{nu,rho}(v, x)
  LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
  for(int nu=0; nu<4; ++nu) {
    LatticePropagator tmp = adj(pl) * gmu5[nu] * pl;
    for(int rho=0; rho<4; ++rho) Gv[nu][rho] = traceS(gL[rho] * tmp) * exp_factor; 
  }

  vector<vector<LatticeColourMatrix>> Cv = conv_with_E_typeI(Fu, env.M_K, max_uv_sep);

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int nu=0; nu<4; ++nu) {
    for(int rho=0; rho<4; ++rho) {
      rst_Q1 += trace(Gv[nu][rho]) * trace(Cv[nu][rho]);
      rst_Q2 += trace(Gv[nu][rho] * Cv[nu][rho]);
    }
  }
}


void typeI_D2a(const std::vector<int> &x, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                     Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                     const LatticePropagator &pl, const LatticePropagator &ps, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagatorSite wl_x; // L(x, tK)
  peekSite(wl_x, wl[tK], x);

  vector<vector<LatticeColourMatrix>> Fu(4, vector<LatticeColourMatrix>(4, env.grid)); // F_{mu,rho}(u, x)
  for(int mu=0; mu<4; ++mu) {
    LatticePropagator A = wl_x * adj(ws[tK]) * gmu5[mu] * ps;
    A = A - adj(A);              //  add contribution from K0_bar
    for(int rho=0; rho<4; ++rho) Fu[mu][rho] = traceS(gL[rho] * A); 
  }

  vector<vector<LatticeColourMatrix>> Gv(4, vector<LatticeColourMatrix>(4, env.grid)); // G_{nu,rho}(v, x)
  LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
  for(int nu=0; nu<4; ++nu) {
    LatticePropagator tmp = adj(pl) * gmu5[nu] * pl;
    for(int rho=0; rho<4; ++rho) Gv[nu][rho] = traceS(gL[rho] * tmp) * exp_factor; 
  }

  vector<vector<LatticeColourMatrix>> Cv = conv_with_E_typeI(Fu, env.M_K, max_uv_sep);

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int nu=0; nu<4; ++nu) {
    for(int rho=0; rho<4; ++rho) {
      rst_Q1 += trace(Gv[nu][rho]) * trace(Cv[nu][rho]);
      rst_Q2 += trace(Gv[nu][rho] * Cv[nu][rho]);
    }
  }
}


void typeI_D1b(const std::vector<int> &x, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                     Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                     const LatticePropagator &pl, const LatticePropagator &ps, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagatorSite ws_x; // H(x, tK)
  peekSite(ws_x, ws[tK], x);

  vector<vector<LatticeColourMatrix>> Fv(4, vector<LatticeColourMatrix>(4, env.grid)); // F_{mu,rho}(u, x)
  LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
  for(int nu=0; nu<4; ++nu) {
    LatticePropagator A = adj(pl) * gmu5[nu] * wl[tK] * adj(ws_x);
    A = A - adj(A);              //  add contribution from K0_bar
    for(int rho=0; rho<4; ++rho) Fv[nu][rho] = traceS(gL[rho] * A) * exp_factor; 
  }

  vector<vector<LatticeColourMatrix>> Gu(4, vector<LatticeColourMatrix>(4, env.grid)); // G_{nu,rho}(v, x)
  for(int mu=0; mu<4; ++mu) {
    LatticePropagator tmp = adj(pl) * gmu5[mu] * pl;
    for(int rho=0; rho<4; ++rho) Gu[mu][rho] = traceS(gL[rho] * tmp); 
  }

  vector<vector<LatticeColourMatrix>> Cv = conv_with_E_typeI(Gu, env.M_K, max_uv_sep);

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int nu=0; nu<4; ++nu) {
    for(int rho=0; rho<4; ++rho) {
      rst_Q1 += trace(Cv[nu][rho]) * trace(Fv[nu][rho]) ;
      rst_Q2 += trace(Cv[nu][rho] * Fv[nu][rho]);
    }
  }
}


void typeI_D2b(const std::vector<int> &x, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                     Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                     const LatticePropagator &pl, const LatticePropagator &ps, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagatorSite wl_x; // L(x, tK)
  peekSite(wl_x, wl[tK], x);

  vector<vector<LatticeColourMatrix>> Fv(4, vector<LatticeColourMatrix>(4, env.grid)); // F_{mu,rho}(u, x)
  LatticeComplex exp_factor = exp_v0_tK(env.grid, tK, env.M_K);
  for(int nu=0; nu<4; ++nu) {
    LatticePropagator A = wl_x * adj(ws[tK]) * gmu5[nu] * ps;
    A = A - adj(A);              //  add contribution from K0_bar
    for(int rho=0; rho<4; ++rho) Fv[nu][rho] = traceS(gL[rho] * A) * exp_factor; 
  }

  vector<vector<LatticeColourMatrix>> Gu(4, vector<LatticeColourMatrix>(4, env.grid)); // G_{nu,rho}(v, x)
  for(int mu=0; mu<4; ++mu) {
    LatticePropagator tmp = adj(pl) * gmu5[mu] * pl;
    for(int rho=0; rho<4; ++rho) Gu[mu][rho] = traceS(gL[rho] * tmp); 
  }

  vector<vector<LatticeColourMatrix>> Cv = conv_with_E_typeI(Gu, env.M_K, max_uv_sep);

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int nu=0; nu<4; ++nu) {
    for(int rho=0; rho<4; ++rho) {
      rst_Q1 += trace(Cv[nu][rho]) * trace(Fv[nu][rho]) ;
      rst_Q2 += trace(Cv[nu][rho] * Fv[nu][rho]);
    }
  }
}


}
