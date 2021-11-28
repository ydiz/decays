#pragma once

#include "kaon/kaon.h"
#include "../control_R.h"

namespace Grid {

using namespace std;

template <class T>
Lattice<T> my_3D_forward_FFT_3D(const Lattice<T> &lat) {

  GridBase *grid = lat.Grid();
  Lattice<T> lat_fft(grid);
  FFT theFFT((GridCartesian *)grid);
  Coordinate mask({1,1,1,0});  // 3D FFT
  theFFT.FFT_dim_mask(lat_fft, lat, mask, FFT::forward);
  return lat_fft;
}

template <class T>
Lattice<T> my_3D_backward_FFT_3D(const Lattice<T> &lat) {

  GridBase *grid = lat.Grid();
  Lattice<T> lat_fft(grid);
  FFT theFFT((GridCartesian *)grid);
  Coordinate mask({1,1,1,0});  // 3D FFT
  theFFT.FFT_dim_mask(lat_fft, lat, mask, FFT::backward);
  return lat_fft;
}

// generate table3d[xt][|vec{x}|][vt]; initialize it to 0
vector<vector<vector<Complex>>> control_R_initialize_table3d(const Coordinate &fdims) { 

    const int T = fdims[3];
    int X = fdims[0], Y = fdims[1], Z = fdims[2];
    int num_R = 1 + int(calc_3d_vec_norm(X/2, Y/2, Z/2));

    vector<vector<vector<Complex>>> table3d(T);  
    for(int xt=0; xt<T; ++xt) {
      table3d[xt].resize(num_R);
      for(int R=0; R<num_R; ++R) {
        table3d[xt][R].resize(T, 0.);
      }
    }
    return table3d;
}


// table3d[xt][R][vt] = sum_{vec{v}} trace( f(xt, vec{x}=vec{v}+R) * g(vt, vec{v}) )
vector<vector<vector<Complex>>> product_and_sum_by_xt_R_vt(const LatticePropagator &fx, const LatticePropagator &gv) {
  LatticePropagator fx_fft = my_3D_forward_FFT_3D(fx);  // Fourier[f]
  LatticePropagator gv_conj = conjugate(gv);
  LatticePropagator gv_fft = conjugate(my_3D_forward_FFT_3D(gv_conj));  // Fourier[f^*]^*

  vector<vector<vector<Complex>>> table3d = control_R_initialize_table3d(fx.Grid()->_fdimensions); // table3d[xt][|vec{x}|][vt] // resize and set initial values to 0

  int T = table3d.size();
  for(int delta_t=0; delta_t<table3d.size(); ++delta_t) { // delta_t :  vt = xt + delta_t

    if(delta_t>=5 && delta_t<=T-4) continue; // only calculate v_t - xt in [-3, 4] // IMPORTANT: to improve performance

    LatticePropagator tmp = fx_fft * Cshift(gv_fft, Tdir, delta_t); // f(x) * g(v = x+delta_t)
    tmp = my_3D_backward_FFT_3D(tmp);
    LatticeComplex amplitude = trace(tmp);

    int num_R = table3d[0].size();
    vector<vector<Complex>> amplitude_by_xt_R = sum_T_R(amplitude, {0,0,0,0});
    for(int xt=0; xt<T; ++xt) {
      int vt = (xt + delta_t) % T;
      for(int R=0; R<num_R; ++R) {
        table3d[xt][R][vt] = amplitude_by_xt_R[xt][R];
      }
    }
  }

  return table3d;
}



// EM factor on the left half lattice; if r_t > 0 or r_t < -max_uv_sep, Euv = 0; If r_t == 0, do not change Euv; If -max_uv_sep<=r_t<0, multiply Euv by 2.
// The factor 2 / M_K^4 is included.
void EM_factor_left_half(LatticeKGG &lat, const vector<int> &u, double M_K, int max_uv_sep) { // max sep is maximum allowed |u0-v0|

  lat = Zero();

  const int T = lat.Grid()->_fdimensions[3];

  autoView(lat_v, lat, CpuWrite);
  parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++) {

    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat.Grid()->_fdimensions);  // smod: if x > T/2, map x to x-T

    if(gcoor[Tdir] > 0 || gcoor[Tdir] < -max_uv_sep) continue; // Non-zero only for - max_uv_sep <= t <= 0
    else {  

      double w = sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2]);

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

  // cout << lat << endl;
  // exit(0);
}

vector<vector<vector<Complex>>> operator+(const vector<vector<vector<Complex>>> &v1, const vector<vector<vector<Complex>>> &v2) {

  vector<vector<vector<Complex>>> rst(v1.size());
  for(int i=0; i<v1.size(); ++i) {
    rst[i].resize(v1[i].size());
    for(int j=0; j<v1[i].size(); ++j) {
      rst[i][j].resize(v1[i][j].size());
      for(int k=0; k<v1[i][j].size(); ++k) {
        rst[i][j][k] = v1[i][j][k] + v2[i][j][k];
      }
    }
  }
  return rst;
}

// Not calculating sBar_d_T2D1b_xt_R_vt for now
void typeIV_D1b(const vector<int> &u, int tK, vector<vector<vector<Complex>>> &rst_Q1_xt_R_vt, 
              vector<vector<vector<Complex>>> &rst_Q2_xt_R_vt, vector<vector<vector<Complex>>> &sBar_d_T2D1b_xt_R_vt,
              Env &env, const vector<LatticePropagator> &wl, const vector<LatticePropagator> &ws, 
              const LatticePropagator &pl, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];
  Coordinate fdims = env.grid->_fdimensions;

  rst_Q1_xt_R_vt = control_R_initialize_table3d(fdims);
  rst_Q2_xt_R_vt = control_R_initialize_table3d(fdims);
  sBar_d_T2D1b_xt_R_vt = control_R_initialize_table3d(fdims);

  int num_R = rst_Q1_xt_R_vt[0].size();

  LatticeKGG Euv(env.grid);
  EM_factor_left_half(Euv, u, env.M_K, max_uv_sep); // EM factor centered at u

  // g = - \sum_v gmu g5 L(v, u)^dagger gnu g5 L(v, tK) E_munu(r=u-v)
  LatticePropagator g_v(env.grid);  g_v = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      g_v += gmu5[mu] * adj(pl) * gmu5[nu] * wl[tK] * Euv_munu;
    }
  }
  g_v = - g_v; // g has a minus sign

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

  rst_Q1_xt_R_vt = product_and_sum_by_xt_R_vt(f_Q1_K, g_v) + product_and_sum_by_xt_R_vt(f_Q1_Kbar, adj(g_v)); // trace(f_Q1_K_xt[xt] * g_v[vt] + f_Q1_Kbar_xt[xt] * adj(g_v[vt]))
  rst_Q2_xt_R_vt = product_and_sum_by_xt_R_vt(f_Q2_K, g_v) + product_and_sum_by_xt_R_vt(f_Q2_Kbar, adj(g_v)); // trace(f_Q2_K_xt[xt] * g_v[vt] + f_Q2_Kbar_xt[xt] * adj(g_v[vt]))
  // sBar_d_T2D1b_xt_R_vt = product_and_sum_by_xt_R_vt(f_sBar_d, g_v);
  
  int tsep = (u[3] - tK + T) % T; 
  double exp_factor = std::exp(env.M_K * tsep);
  for(int xt=0; xt<T; ++xt) {
    for(int R=0; R<num_R; ++R) {
      for(int vt=0; vt<T; ++vt) {
        // sBar_d_T2D1b_xt_R_vt[xt][R][vt] = 2. * real(sBar_d_T2D1b_xt_R_vt[xt][R][vt]); // add contribution of K0bar
        // sBar_d_T2D1b_xt_R_vt[xt][R][vt] *= exp_factor;

        rst_Q1_xt_R_vt[xt][R][vt] *= exp_factor;  // multiply all amplitudes by exp(M_K * (u0 - tK))
        rst_Q2_xt_R_vt[xt][R][vt] *= exp_factor;
      }
    }
  }

}



void typeIV_D2b(const vector<int> &u, int tK, vector<vector<vector<Complex>>> &rst_Q1_xt_R_vt, 
              vector<vector<vector<Complex>>> &rst_Q2_xt_R_vt, vector<vector<vector<Complex>>> &sBar_d_T2D2b_xt_R_vt,
              Env &env, const vector<LatticePropagator> &wl, const vector<LatticePropagator> &ws, 
              const LatticePropagator &ps, const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];
  Coordinate fdims = env.grid->_fdimensions;

  rst_Q1_xt_R_vt = control_R_initialize_table3d(fdims);
  rst_Q2_xt_R_vt = control_R_initialize_table3d(fdims);
  sBar_d_T2D2b_xt_R_vt = control_R_initialize_table3d(fdims);

  int num_R = rst_Q1_xt_R_vt[0].size();

  LatticeKGG Euv(env.grid);
  EM_factor_left_half(Euv, u, env.M_K, max_uv_sep); // EM factor centered at u

  // g = - \sum_v gmu g5 L(v, u)^dagger gnu g5 L(v, tK) E_munu(r=u-v)
  LatticePropagator g_v(env.grid);  g_v = Zero();
  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // For these indicies, Euv_munu = 0
      LatticeComplex Euv_munu = PeekIndex<LorentzIndex>(Euv, mu, nu);
      g_v += adj(ws[tK]) * gmu5[nu] * ps * gmu5[mu] * Euv_munu;
    }
  }
  g_v = - g_v; // g has a minus sign

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

  rst_Q1_xt_R_vt = product_and_sum_by_xt_R_vt(f_Q1_K, g_v) + product_and_sum_by_xt_R_vt(f_Q1_Kbar, adj(g_v)); // trace(f_Q1_K_xt[xt] * g_v[vt] + f_Q1_Kbar_xt[xt] * adj(g_v[vt]))
  rst_Q2_xt_R_vt = product_and_sum_by_xt_R_vt(f_Q2_K, g_v) + product_and_sum_by_xt_R_vt(f_Q2_Kbar, adj(g_v)); // trace(f_Q2_K_xt[xt] * g_v[vt] + f_Q2_Kbar_xt[xt] * adj(g_v[vt]))
  // sBar_d_T2D2b_xt_R_vt = product_and_sum_by_xt_R_vt(f_sBar_d, g_v);

  int tsep = (u[3] - tK + T) % T;
  double exp_factor = std::exp(env.M_K * tsep);
  for(int xt=0; xt<T; ++xt) {
    for(int R=0; R<num_R; ++R) {
      for(int vt=0; vt<T; ++vt) {
        // sBar_d_T2D2b_xt_R_vt[xt][R][vt] = 2. * real(sBar_d_T2D2b_xt_R_vt[xt][R][vt]); // add contribution of K0bar
        // sBar_d_T2D2b_xt_R_vt[xt][R][vt] *= exp_factor;

        rst_Q1_xt_R_vt[xt][R][vt] *= exp_factor;  // multiply all amplitudes by exp(M_K * (u0 - tK))
        rst_Q2_xt_R_vt[xt][R][vt] *= exp_factor;
      }
    }
  }

}





}

