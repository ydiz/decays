
// the result is the form factor except for the coefficient C1, C2 and the coefficient of each individual diagram

#include "../kaon.h"
#include "amplitude/form_factor.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;





LatticeKGG calc_leptonic_with_coef(double M_K, const std::vector<int> &v, GridCartesian *grid) {  // return L_munu(u, v), where v is a fixed point.

  LatticeKGG lep(grid);

  form_factor_integrand(lep, M_K);

  double lep_coef = 2. / std::pow(M_K, 4);
  // double hadron_coef = env.Z_V * env.Z_V * 2. * env.M_K / env.N_K;  // Note: I did not multiply hadronic coefficient
  
  lep = lep * lep_coef;

  for(int mu=0; mu<4; ++mu) lep = Cshift(lep, mu, -v[mu]); // shift v to origin // FIXME: check it is right to shift -v[mu]

  return lep;
}


void restrictTimeRange(LatticePropagator &lat, int vt)  { // The allowed interval of u is: vt <= ut <= vt+16

  const int T = lat.Grid()->_fdimensions[3];

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    int ut = gcoor[3];

    int dist;   // dist = u_t - v_t, taking periodic bundary condition into account
    if(abs(ut - vt)<=T/2) dist = ut - vt;
    else if(abs(ut + T - vt) <= T/2) dist = ut + T - vt;
    else dist = ut - T - vt;

    LatticePropagatorSite m;
    if(dist==0) {}                      // if u_t == v_t , do nothing
    else if(dist > 0 && dist <= T/4) {  // u_t>v_t && u_t - v_t <= T/4  // if u is on the right of v, multiple it by 2
      peekLocalSite(m, lat, lcoor);
      m = 2. * m;
      pokeLocalSite(m, lat, lcoor);
    }   
    else {                           // else, set to 0
      m = Zero(); 
      pokeLocalSite(m, lat, lcoor);
    }   

  }); 
}




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
    LatticePropagator Lxx = env.get_Lxx();

    int num_timeSlices = T/2 - tsep3 - tsep2 + 1;
    vector<Complex> rst_Q1_avgSrc(num_timeSlices, 0.), rst_Q2_avgSrc(num_timeSlices, 0.); // average amplitude on each time slice after averaging over 512 point sources

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const auto &v: env.xgs_l) {
      std::cout << "# Point source: " << num_pt_src << std::endl;
      ++num_pt_src;

      int tK = v[3] - tsep;
      if(tK < 0) tK += T;

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v)



      // f1 = \sum_u gnu g5 L(u, v)^dagger gmu g5 L(u, tK)
      LatticePropagator f1_tmp(env.grid);  f1_tmp = Zero();
      LatticeKGG lep = calc_leptonic_with_coef(env.M_K, v, env.grid);
      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          if( mu==3 || nu==3 || mu==nu) continue;  // For these directions, lep_munu = 0
          LatticeComplex lep_munu = PeekIndex<LorentzIndex>(lep, mu, nu);
          f1_tmp += gmu5[nu] * adj(pl) * gmu5[mu] * wl[tK] * lep_munu;
        }
      }
      // // the summand is non-zero only for the region where vt <= ut <= vt+T/4; When ut>vt, multiply by 2
      restrictTimeRange(f1_tmp, v[3]); 
      LatticePropagatorSite f1 = sum(f1_tmp);


      // calculate contraction
      LatticeComplex tmp_Q1(env.grid), tmp_Q2(env.grid);
      tmp_Q1 = Zero(); tmp_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        tmp_Q1 += trace(gL[rho] * Lxx) * ( trace( adj(ws[tK]) * gL[rho] * pl * f1 ) 
                                           - trace( adj(pl) * gL[rho] * ws[tK]  * adj(f1) ) );   // K bar
        tmp_Q2 += trace( adj(ws[tK]) * gL[rho] * Lxx *  gL[rho] * pl * f1 ) 
                  - trace( adj(pl) * gL[rho] * Lxx * gL[rho] * ws[tK] * adj(f1) );
      }

      // sum over each time slice 
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
