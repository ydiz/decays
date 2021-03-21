#pragma once

#include "kaon/kaon.h"

namespace Grid {



void typeIII_D1a(const std::vector<int> &v, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &pl, const LatticePropagator &sequential) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagator f1(env.grid);
  f1 = sequential * g5 * adj(pl);

  LatticePropagator f2(env.grid);
  f2 = wl[tK] * adj(ws[tK]);
  f2 = f2 + adj(f2);             // to incorporate the contribution of K0 bar f2 = wl[tK] * adj(ws[tK]) - ws[tK] * adj(wl[tK]) 

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int rho=0; rho<4; ++rho) {
    rst_Q1 += trace(gL[rho] * f1) *  trace(gL[rho] * f2);
    rst_Q2 += trace(gL[rho] * f1 * gL[rho] * f2);
  }

  int tsep = (v[3] - tK + T) % T;
  rst_Q1 *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
  rst_Q2 *= std::exp(env.M_K * tsep);
}



void typeIII_D1b(const std::vector<int> &v, int tK, LatticeComplex &rst_Q1, LatticeComplex &rst_Q2, 
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &pl, const LatticePropagator &sequential) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagator f1(env.grid);
  f1 = sequential * g5 * adj(pl);
  f1 = adj(f1);               // This line is the only difference between typeIII_D1a and typeIII_D1b

  LatticePropagator f2(env.grid);
  f2 = wl[tK] * adj(ws[tK]);
  f2 = f2 + adj(f2);             // to incorporate the contribution of K0 bar f2 = wl[tK] * adj(ws[tK]) - ws[tK] * adj(wl[tK]) 

  rst_Q1 = Zero(); rst_Q2 = Zero();
  for(int rho=0; rho<4; ++rho) {
    rst_Q1 += trace(gL[rho] * f1) *  trace(gL[rho] * f2);
    rst_Q2 += trace(gL[rho] * f1 * gL[rho] * f2);
  }

  int tsep = (v[3] - tK + T) % T;
  rst_Q1 *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
  rst_Q2 *= std::exp(env.M_K * tsep);
}



}
