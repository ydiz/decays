#pragma once

#include "kaon/kaon.h"

namespace Grid {


// EM factor on the left half lattice; if r_t > 0 or r_t < -max_uv_sep, Euv = 0; If r_t == 0, do not change Euv; If -max_uv_sep<=r_t<0, multiply Euv by 2.
// The factor 2 / M_K^4 is included.
void EM_factor_left_half(LatticeKGG &lat, const std::vector<int> &u, double M_K, int max_uv_sep) { // max sep is maximum allowed |u0-v0|

  lat = Zero();

  const int T = lat.Grid()->_fdimensions[3];

  autoView(lat_v, lat, CpuWrite);
  parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++) {

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat.Grid()->_fdimensions);  // smod: if x > T/2, map x to x-T

    if(gcoor[Tdir] > 0 || gcoor[Tdir] < -max_uv_sep) continue; // Non-zero only for - max_uv_sep <= t <= 0
    else {  

      double w = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2]);

      double val;
      if(w==0) val = 0.;
      else {
        double t = M_K * 0.5 * w;
        double sin_t, cos_t;
        sincos(t, &sin_t, &cos_t);
        val = ( - cos_t * M_K * w + 2 * sin_t) / w / w / w * std::exp(M_K * 0.5 * gcoor[Tdir]);
      }

      if(gcoor[Tdir] < 0) val *= 2.;  // if t < 0, multiply it by 2.

      typename LatticeKGG::vector_object::scalar_object m;
      m = Zero();
      m(0, 1)()() = Complex(val * gcoor[Zdir], 0); 
      m(0, 2)()() = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
      m(1, 2)()() = Complex(val * gcoor[Xdir], 0); 
      m(1, 0)()() = - m(0, 1)()();
      m(2, 0)()() = - m(0, 2)()();
      m(2, 1)()() = - m(1, 2)()();

      pokeLocalSite(m, lat_v, lcoor);
    }
  }

  lat *= (2. / M_K / M_K / M_K / M_K);

  for(int mu=0; mu<4; ++mu) {  // shift center to u
    if(u[mu] != 0) lat = Cshift(lat, mu, -u[mu]);
  }

  // std::cout << lat << std::endl;
  // exit(0);
}




void typeIV_D1b(const std::vector<int> &u, int tK, std::vector<std::vector<Complex>> &rst_Q1_xt_vt, 
                std::vector<std::vector<Complex>> &rst_Q2_xt_vt,
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &pl, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticeKGG Euv(env.grid);
  EM_factor_left_half(Euv, u, env.M_K, max_uv_sep); // EM factor centered at u

  // g = \sum_u gnu g5 L(u, v)^dagger gmu g5 L(u, tK)
  LatticePropagator g_u(env.grid);  g_u = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      g_u += gmu5[mu] * adj(pl) * gmu5[nu] * wl[tK] * Euv_munu;
    }
  }
  g_u = - g_u; // g has a minus sign

  vector<LatticePropagatorSite> g_vt;
  sliceSum(g_u, g_vt, Tdir);
  // LatticePropagatorSite g = sum(g_u);

  // calculate contraction
  LatticePropagator f_Q1_K(env.grid), f_Q1_Kbar(env.grid), f_Q2_K(env.grid), f_Q2_Kbar(env.grid);
  f_Q1_K = Zero(); f_Q1_Kbar = Zero(); f_Q2_K = Zero(); f_Q2_Kbar = Zero();
  for(int rho=0; rho<4; ++rho) {
    f_Q1_K += trace(gL[rho] * Lxx) * adj(ws[tK]) * gL[rho] * pl;
    f_Q1_Kbar += trace(gL[rho] * Lxx) * adj(pl) * gL[rho] * ws[tK];

    f_Q2_K += adj(ws[tK]) * gL[rho] * Lxx *  gL[rho] * pl;
    f_Q2_Kbar += adj(pl) * gL[rho] * Lxx * gL[rho] * ws[tK];
  }

  vector<LatticePropagatorSite> f_Q1_K_xt; sliceSum(f_Q1_K, f_Q1_K_xt, Tdir);
  vector<LatticePropagatorSite> f_Q1_Kbar_xt; sliceSum(f_Q1_Kbar, f_Q1_Kbar_xt, Tdir);
  vector<LatticePropagatorSite> f_Q2_K_xt; sliceSum(f_Q2_K, f_Q2_K_xt, Tdir);
  vector<LatticePropagatorSite> f_Q2_Kbar_xt; sliceSum(f_Q2_Kbar, f_Q2_Kbar_xt, Tdir);

  rst_Q1_xt_vt.resize(T); for(auto &x: rst_Q1_xt_vt) x.resize(T);
  rst_Q2_xt_vt.resize(T); for(auto &x: rst_Q2_xt_vt) x.resize(T);

  int tsep = (u[3] - tK + T) % T;
  for(int xt=0; xt<T; ++xt) {
    for(int vt=0; vt<T; ++vt) {
      rst_Q1_xt_vt[xt][vt] = trace(f_Q1_K_xt[xt] * g_vt[vt] + f_Q1_Kbar_xt[xt] * adj(g_vt[vt]))()()();
      rst_Q2_xt_vt[xt][vt] = trace(f_Q2_K_xt[xt] * g_vt[vt] + f_Q2_Kbar_xt[xt] * adj(g_vt[vt]))()()();

      rst_Q1_xt_vt[xt][vt] *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (u0 - tK))
      rst_Q2_xt_vt[xt][vt] *= std::exp(env.M_K * tsep);
    }
  }

}



void typeIV_D2b(const std::vector<int> &u, int tK, std::vector<std::vector<Complex>> &rst_Q1_xt_vt, 
                std::vector<std::vector<Complex>> &rst_Q2_xt_vt,
                Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                const LatticePropagator &ps, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticeKGG Euv(env.grid);
  EM_factor_left_half(Euv, u, env.M_K, max_uv_sep); // EM factor centered at u

  // g = \sum_u gnu g5 L(u, v)^dagger gmu g5 L(u, tK)
  LatticePropagator g_u(env.grid);  g_u = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      g_u += adj(ws[tK]) * gmu5[nu] * ps * gmu5[mu] * Euv_munu;
    }
  }
  g_u = - g_u; // g has a minus sign

  // LatticePropagatorSite g = sum(g_u);
  vector<LatticePropagatorSite> g_vt;
  sliceSum(g_u, g_vt, Tdir);

  // calculate contraction
  LatticePropagator f_Q1_K(env.grid), f_Q1_Kbar(env.grid), f_Q2_K(env.grid), f_Q2_Kbar(env.grid);
  f_Q1_K = Zero(); f_Q1_Kbar = Zero(); f_Q2_K = Zero(); f_Q2_Kbar = Zero();
  for(int rho=0; rho<4; ++rho) {
    f_Q1_K += trace(gL[rho] * Lxx) * adj(ps) * gL[rho] * wl[tK];
    f_Q1_Kbar += trace(gL[rho] * Lxx) * adj(wl[tK]) * gL[rho] * ps;

    f_Q2_K += adj(ps) * gL[rho] * Lxx *  gL[rho] * wl[tK];
    f_Q2_Kbar += adj(wl[tK]) * gL[rho] * Lxx * gL[rho] * ps;
  }

  vector<LatticePropagatorSite> f_Q1_K_xt; sliceSum(f_Q1_K, f_Q1_K_xt, Tdir);
  vector<LatticePropagatorSite> f_Q1_Kbar_xt; sliceSum(f_Q1_Kbar, f_Q1_Kbar_xt, Tdir);
  vector<LatticePropagatorSite> f_Q2_K_xt; sliceSum(f_Q2_K, f_Q2_K_xt, Tdir);
  vector<LatticePropagatorSite> f_Q2_Kbar_xt; sliceSum(f_Q2_Kbar, f_Q2_Kbar_xt, Tdir);

  rst_Q1_xt_vt.resize(T); for(auto &x: rst_Q1_xt_vt) x.resize(T);
  rst_Q2_xt_vt.resize(T); for(auto &x: rst_Q2_xt_vt) x.resize(T);

  int tsep = (u[3] - tK + T) % T;
  for(int xt=0; xt<T; ++xt) {
    for(int vt=0; vt<T; ++vt) {
      rst_Q1_xt_vt[xt][vt] = trace(f_Q1_K_xt[xt] * g_vt[vt] + f_Q1_Kbar_xt[xt] * adj(g_vt[vt]))()()();
      rst_Q2_xt_vt[xt][vt] = trace(f_Q2_K_xt[xt] * g_vt[vt] + f_Q2_Kbar_xt[xt] * adj(g_vt[vt]))()()();

      rst_Q1_xt_vt[xt][vt] *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (u0 - tK))
      rst_Q2_xt_vt[xt][vt] *= std::exp(env.M_K * tsep);
    }
  }


  // rst_Q1 = trace(f_Q1_K * g - f_Q1_Kbar * adj(g));
  // rst_Q2 = trace(f_Q2_K * g - f_Q2_Kbar * adj(g));
  //
  // int tsep = (v[3] - tK + T) % T;
  // rst_Q1 *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
  // rst_Q2 *= std::exp(env.M_K * tsep);
}





}

