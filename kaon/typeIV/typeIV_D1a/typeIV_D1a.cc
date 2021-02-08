#include "kaon/kaon.h"


namespace Grid {

// For each v and tK, return amplitude as a function of x_t
void typeIV_diagram_1a(const std::vector<int> &v, int tK, std::vector<Complex> &rst_Q1_xt, std::vector<Complex> &rst_Q2_xt, 
                     Env &env, const std::vector<LatticePropagator> &wl, const std::vector<LatticePropagator> &ws, 
                     const LatticePropagator &Lxx, int max_uv_sep) {
  using namespace std;
  const int T = env.grid->_fdimensions[3];

  LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) 

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
  LatticePropagator f_Q1_K(env.grid), f_Q1_Kbar(env.grid), f_Q2_K(env.grid), f_Q2_Kbar(env.grid);
  f_Q1_K = Zero(); f_Q1_Kbar = Zero(); f_Q2_K = Zero(); f_Q2_Kbar = Zero();
  for(int rho=0; rho<4; ++rho) {
    f_Q1_K += trace(gL[rho] * Lxx) * adj(ws[tK]) * gL[rho] * pl;
    f_Q1_Kbar += trace(gL[rho] * Lxx) * adj(pl) * gL[rho] * ws[tK];

    f_Q2_K += adj(ws[tK]) * gL[rho] * Lxx *  gL[rho] * pl;
    f_Q2_Kbar += adj(pl) * gL[rho] * Lxx * gL[rho] * ws[tK];
  }

  vector<LatticePropagatorSite> f_Q1_K_xt, f_Q1_Kbar_xt, f_Q2_K_xt, f_Q2_Kbar_xt;
  sliceSum(f_Q1_K, f_Q1_K_xt, Tdir);
  sliceSum(f_Q1_Kbar, f_Q1_Kbar_xt, Tdir);
  sliceSum(f_Q2_K, f_Q2_K_xt, Tdir);
  sliceSum(f_Q2_Kbar, f_Q2_Kbar_xt, Tdir);

  rst_Q1_xt.resize(T); rst_Q2_xt.resize(T);
  for(int i=0; i<T; ++i) {
    rst_Q1_xt[i] = trace(f_Q1_K_xt[i] * g - f_Q1_Kbar_xt[i] * adj(g))()()();
    rst_Q2_xt[i] = trace(f_Q2_K_xt[i] * g - f_Q2_Kbar_xt[i] * adj(g))()()();
  }

  int tsep = (v[3] - tK + T) % T;
  for(int i=0; i<T; ++i) {
    rst_Q1_xt[i] *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
    rst_Q2_xt[i] *= std::exp(env.M_K * tsep);
  }
}

}

using namespace std;
using namespace Grid;

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  // zyd_init_Grid_Qlattice(argc, argv);

  int target_traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, target_traj);
  }
  else {
    std::cout << "traj not specified; exiting" << std::endl;
    assert(0);
  }
  int traj_start = target_traj;
  int traj_end = target_traj;
  int traj_sep = 10;

  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

#ifdef CUTH_FREE_FIELD
  int tsep = 1;
  int max_uv_sep = 3;
  Env env("FreeField_8nt8");
  // std::vector<int> allowed_x0_tK = {1,2,3}; // values of x0-tK // FIXME: may change this later.
#else
  int tsep = 12;  // tv = tK + tsep
  int max_uv_sep = 16;
  // std::vector<int> allowed_x0_tK = {4,6,8,10,12,14}; // values of x0-tK // FIXME: may change this later.
  Env env("24ID");
#endif

  env.N_pt_src = -1;  

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');
    LatticePropagator Lxx = env.get_Lxx();

    vector<Complex> rst_Q1_xt_avgSrc(T, 0.), rst_Q2_xt_avgSrc(T, 0.); // average amplitude on each time slice of x after averaging over 512 point sources // Note: index 0 is the position of point v

    vector<vector<vector<Complex>>> table3d; // A[tK, x0, v0]
    vector<int> n_psrc_time_slice(T, 0); // number of point source on each time slice

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const vector<int> &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      std::cout << "Point source: " << v << std::endl;
      ++num_pt_src;
      n_psrc_time_slice[v[3]] += 1;

      int tK = (v[3] - tsep + T) % T;

      vector<Complex> rst_Q1_xt(T), rst_Q2_xt(T);
      typeIV_diagram_1a(v, tK, rst_Q1_xt, rst_Q2_xt, env, wl, ws, Lxx, max_uv_sep);

      for(int i=0; i<T; ++i) {
        rst_Q1_xt_avgSrc[i] += rst_Q1_xt[(i+v[3])%T]; // rst_Q1_avgSrc[0] is the time slice of point v
        rst_Q2_xt_avgSrc[i] += rst_Q2_xt[(i+v[3])%T];
      }

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;


    for(int i=0; i<T; ++i) {
      rst_Q1_xt_avgSrc[i] /= double(num_pt_src);
      rst_Q2_xt_avgSrc[i] /= double(num_pt_src);
    }

    std::cout << "traj [" << traj << "] amplitude Q1 by time slice: " << rst_Q1_xt_avgSrc << std::endl;
    std::cout << "traj [" << traj << "] amplitude Q2 by time slice: " << rst_Q2_xt_avgSrc << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
