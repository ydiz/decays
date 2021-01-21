
// the result is the form factor except for the coefficient C1, C2 and the coefficient of each individual diagram

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

  double hadron_coef = env.Z_V * env.Z_V * 2. * env.M_K / env.N_K;  // coefficient of hadronic part

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    int num_timeSlices = T/2 - tsep3 - tsep2 + 1;
    vector<Complex> rst_Q1_avgSrc(num_timeSlices, 0.), rst_Q2_avgSrc(num_timeSlices, 0.); // average amplitude on each time slice after averaging over 512 point sources

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const auto &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      ++num_pt_src;

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) 

      LatticePropagator sequential = env.get_sequential(v);

      int tK = v[3] - tsep;
      if(tK < 0) tK += T;

      LatticePropagator f1(env.grid);
      f1 = sequential * g5 * adj(pl);

      LatticePropagator f2(env.grid);
      f2 = wl[tK] * adj(ws[tK]);
      f2 = f2 - adj(f2);             // to incorporate the contribution of K0 bar  f2 = wl[tK] * adj(ws[tK]) - ws[tK] * adj(wl[tK]) 

      LatticeComplex tmp_Q1(env.grid), tmp_Q2(env.grid);
      tmp_Q1 = Zero(); tmp_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        tmp_Q1 += trace(gL[rho] * f1) *  trace(gL[rho] * f2);
        tmp_Q2 += trace(gL[rho] * f1 * gL[rho] * f2);
      }

      int lower_bound = tK + tsep2, upper_bound = tK + T/2 - tsep3; // both lower_bound and upper_bound can be greater than T
      Sum_Interval_TimeSlice sum_interval(lower_bound, upper_bound, T); 
      vector<LatticeComplexSite> rst_Q1_tmp = sum_interval(tmp_Q1);
      vector<LatticeComplexSite> rst_Q2_tmp = sum_interval(tmp_Q2);
      vector<Complex> rst_Q1; for(LatticeComplexSite x: rst_Q1_tmp) rst_Q1.push_back(x()()());
      vector<Complex> rst_Q2; for(LatticeComplexSite x: rst_Q2_tmp) rst_Q2.push_back(x()()());
      assert(rst_Q1.size() == rst_Q1_avgSrc.size());
      assert(rst_Q2.size() == rst_Q2_avgSrc.size());

      for(int i=0; i<rst_Q1.size(); ++i) rst_Q1_avgSrc[i] += rst_Q1[i];
      for(int i=0; i<rst_Q2.size(); ++i) rst_Q2_avgSrc[i] += rst_Q2[i];


    } // end of point source loop

    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;


    for(int i=0; i<rst_Q1_avgSrc.size(); ++i) {
      rst_Q1_avgSrc[i] *= std::exp(env.M_K * tsep);
      rst_Q2_avgSrc[i] *= std::exp(env.M_K * tsep);

      rst_Q1_avgSrc[i] /= double(num_pt_src);
      rst_Q2_avgSrc[i] /= double(num_pt_src);

      rst_Q1_avgSrc[i] *= hadron_coef;
      rst_Q2_avgSrc[i] *= hadron_coef;
    }

    std::cout << "traj [" << traj << "] amplitude Q1 by time slice: " << rst_Q1_avgSrc << std::endl;
    std::cout << "traj [" << traj << "] amplitude Q2 by time slice: " << rst_Q2_avgSrc << std::endl;

    Complex total_amplitude_Q1 = 0., total_amplitude_Q2 = 0.;
    for(Complex x: rst_Q1_avgSrc) total_amplitude_Q1 += x;
    for(Complex x: rst_Q2_avgSrc) total_amplitude_Q2 += x;

    std::cout << "traj [" << traj << "] total amplitude Q1: " << total_amplitude_Q1 << std::endl;
    std::cout << "traj [" << traj << "] total amplitude Q2: " << total_amplitude_Q2 << std::endl;

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
