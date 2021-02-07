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
  Env env("FreeField_8nt8");
#else
  int tsep = 12;  // tv = tK + tsep
  Env env("24ID");
#endif

  env.N_pt_src = -1;  

  // double hadron_coef = env.Z_V * env.Z_V * 2. * env.M_K / env.N_K;  // coefficient of hadronic part

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    // int num_timeSlices = T/2 - tsep3 - tsep2 + 1;
    // vector<Complex> rst_Q1_avgSrc(num_timeSlices, 0.), rst_Q2_avgSrc(num_timeSlices, 0.); // average amplitude on each time slice after averaging over 512 point sources

    vector<Complex> rst_Q1_xt_avgSrc(T, 0.), rst_Q2_xt_avgSrc(T, 0.); // average amplitude on each time slice after averaging over 512 point sources // Note: index 0 is the position of point v

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const auto &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      ++num_pt_src;

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) 

      LatticePropagator sequential = env.get_sequential(v);

      int tK = (v[3] - tsep + T) % T;

      LatticePropagator f1(env.grid);
      f1 = sequential * g5 * adj(pl);

      LatticePropagator f2(env.grid);
      f2 = wl[tK] * adj(ws[tK]);
      f2 = f2 - adj(f2);             // to incorporate the contribution of K0 bar f2 = wl[tK] * adj(ws[tK]) - ws[tK] * adj(wl[tK]) 

      // {  // for test
      //   LatticeComplex tmp(env.grid);
      //   // tmp = trace(g5 * f1);
      //   tmp = trace(f1);
      //   std::cout << tmp << std::endl;
      //   exit(0);
      // }

      LatticeComplex rst_Q1(env.grid), rst_Q2(env.grid);
      rst_Q1 = Zero(); rst_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        rst_Q1 += trace(gL[rho] * f1) *  trace(gL[rho] * f2);
        // rst_Q1 += trace(g5 * f1) *  trace(gL[rho] * f2);
        rst_Q2 += trace(gL[rho] * f1 * gL[rho] * f2);
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
