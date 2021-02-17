#include "kaon/kaon.h"

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
#else
  int tsep = 12;  // tv = tK + tsep
  int max_uv_sep = 16;
  Env env("24ID");
#endif

  env.N_pt_src = -1;  

  // double hadron_coef = env.Z_V * env.Z_V * 2. * env.M_K / env.N_K;  // coefficient of hadronic part

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');
    LatticePropagator Lxx = env.get_Lxx();

    // vector<Complex> rst_K0_Q1_xt_avgSrc(T, 0.), rst_K0_Q2_xt_avgSrc(T, 0.), rst_K0bar_Q1_xt_avgSrc(T, 0.), rst_K0bar_Q2_xt_avgSrc(T, 0.); // average amplitude on each time slice of x after averaging over 512 point sources // Note: index 0 is the position of point v
    vector<Complex> rst_Q1_xt_avgSrc(T, 0.), rst_Q2_xt_avgSrc(T, 0.); // average amplitude on each time slice of x after averaging over 512 point sources // Note: index 0 is the position of point v

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const auto &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      std::cout << "Point source: " << v << std::endl;
      ++num_pt_src;

      LatticePropagator ps = env.get_point(v, 's'); // ps = H(x, v) 

      int tK = (v[3] - tsep + T) % T;

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
      LatticeComplex rst_Q1(env.grid), rst_Q2(env.grid);
      rst_Q1 = Zero(); rst_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        rst_Q1 += trace(gL[rho] * Lxx) * ( trace( adj(ps) * gL[rho] * wl[tK] * g ) 
                                           - trace( adj(wl[tK]) * gL[rho] * ps * adj(g) ) );   // K bar
        rst_Q2 += trace( adj(ps) * gL[rho] * Lxx *  gL[rho] * wl[tK] * g) 
                  - trace( adj(wl[tK]) * gL[rho] * Lxx * gL[rho] * ps * adj(g) );
      }


      vector<LatticeComplexSite> rst_Q1_xt, rst_Q2_xt;
      sliceSum(rst_Q1, rst_Q1_xt, Tdir);
      sliceSum(rst_Q2, rst_Q2_xt, Tdir);
      for(int i=0; i<T; ++i) {
        rst_Q1_xt_avgSrc[i] += rst_Q1_xt[(i+v[3])%T]()()(); // rst_Q1_avgSrc[0] is the time slice of point v
        rst_Q2_xt_avgSrc[i] += rst_Q2_xt[(i+v[3])%T]()()();
      }

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;


    for(int i=0; i<T; ++i) {
      rst_Q1_xt_avgSrc[i] *= std::exp(env.M_K * tsep);  // multiply all amplitudes by exp(M_K * (v0 - tK))
      rst_Q2_xt_avgSrc[i] *= std::exp(env.M_K * tsep);

      rst_Q1_xt_avgSrc[i] /= double(num_pt_src);
      rst_Q2_xt_avgSrc[i] /= double(num_pt_src);

      // rst_Q1_avgSrc[i] *= hadron_coef; // Do not multiply hadron_coef in C++
      // rst_Q2_avgSrc[i] *= hadron_coef;
    }

    std::cout << "traj [" << traj << "] amplitude Q1 by time slice: " << rst_Q1_xt_avgSrc << std::endl;
    std::cout << "traj [" << traj << "] amplitude Q2 by time slice: " << rst_Q2_xt_avgSrc << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}