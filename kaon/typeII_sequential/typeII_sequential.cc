
#include "../kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char* argv[])
{
  // Grid_init(&argc, &argv);
  zyd_init_Grid_Qlattice(argc, argv);

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

  // change those parameters
  int tsep = 12;  // tv = tK + tsep
  int tsep2 = 6;  // tx >= tK + tsep2
  int tsep3 = 4;  // tx <= tK + T/2 - tsep3

  Env env("24ID");
  env.N_pt_src = -1;  

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    Complex amplitude_Q1_allsrc = 0., amplitude_Q2_allsrc = 0.;

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const auto &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      ++num_pt_src;

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) 

      LatticePropagator sequential = env.get_sequential(v, "typeII");

      int tK = v[3] - tsep;
      if(tK < 0) tK += T;

      LatticePropagator f1(env.grid), f2(env.grid);
      f1 = sequential * g5 * adj(pl);
      f2 = wl[tK] * adj(ws[tK]);

      LatticeComplex tmp_Q1(env.grid), tmp_Q2(env.grid);
      tmp_Q1 = Zero(); tmp_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        tmp_Q1 += trace(gL[rho] * f1) *  trace(gL[rho] * f2);
        tmp_Q2 += trace(gL[rho] * f1 * gL[rho] * f2);
      }

      int lower_bound = tK + tsep2, upper_bound = tK + T/2 - tsep3; // both lower_bound and upper_bound can be greater than T
      Sum_Interval sum_interval(lower_bound, upper_bound, T);   

      Complex rst_Q1 = sum_interval(tmp_Q1)()()();
      Complex rst_Q2 = sum_interval(tmp_Q2)()()();

      rst_Q1 *= std::exp(env.M_K * tsep);
      rst_Q2 *= std::exp(env.M_K * tsep);

      amplitude_Q1_allsrc += rst_Q1;
      amplitude_Q2_allsrc += rst_Q2;

    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;

    amplitude_Q1_allsrc /= double(num_pt_src);
    amplitude_Q2_allsrc /= double(num_pt_src);
    std::cout << "amplitude Q1: " << amplitude_Q1_allsrc << std::endl;
    std::cout << "amplitude Q2: " << amplitude_Q2_allsrc << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
