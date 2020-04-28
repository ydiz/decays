#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


struct Sum_Interval {
  int lower_bound, upper_bound, T;

  Sum_Interval(int _lower_bound, int _upper_bound, int _T) : lower_bound(_lower_bound), upper_bound(_upper_bound), T(_T) {
    assert(lower_bound >=0 && upper_bound>=0 && upper_bound >= lower_bound); // upper_bound can be greater than T
  };

  LatticePropagatorSite operator()(const LatticePropagator &lat) const { 
    std::vector<LatticePropagatorSite> lat_sliceSum;
    sliceSum(lat, lat_sliceSum, Tdir);

    LatticePropagatorSite rst = Zero();
    for(int t=lower_bound; t<=upper_bound; ++t) {
      rst += (t>=T ? lat_sliceSum[t % T] : lat_sliceSum[t]);
    }
    return rst;
  }
};



std::vector<int> gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  // FIXME: change those parameters
  int tsep = 10;
  int tsep2 = 6;
  int tsep3 = 6;

  Env env(gcoor, "24ID");
  // init_para(argc, argv, env);
  env.N_pt_src = 1;  // FIXME: keep only one point

  env.grid->show_decomposition();

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    if(env.N_pt_src != -1) env.xgs_s.erase(env.xgs_s.begin() + env.N_pt_src, env.xgs_s.end());

    LatticePropagator Lxx(env.grid);
    Lxx = 1.0; // FIXME: Should use A2A propagator! I am setting one for test.

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    // save average of all point sources
    LatticePGG rst_D1Q1K_allsrc(env.grid), rst_D1Q2K_allsrc(env.grid), rst_D1Q1Kbar_allsrc(env.grid), rst_D1Q2Kbar_allsrc(env.grid); 
    LatticePGG rst_D2Q1K_allsrc(env.grid), rst_D2Q2K_allsrc(env.grid), rst_D2Q1Kbar_allsrc(env.grid), rst_D2Q2Kbar_allsrc(env.grid); 
    LatticePGG rst_sBar_d_D1_allsrc(env.grid), rst_sBar_d_D2_allsrc(env.grid); 
    vector<LatticePGG*> rst_vec_allsrc= {&rst_D1Q1K_allsrc, &rst_D1Q2K_allsrc, &rst_D1Q1Kbar_allsrc, &rst_D1Q2Kbar_allsrc, 
                                   &rst_D2Q1K_allsrc, &rst_D2Q2K_allsrc, &rst_D2Q1Kbar_allsrc, &rst_D2Q2Kbar_allsrc, 
                                   &rst_sBar_d_D1_allsrc, &rst_sBar_d_D2_allsrc};
    for(auto rst: rst_vec_allsrc) *rst = Zero();

    for(const auto &v: env.xgs_s) {

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) or L(u, v)
      LatticePropagator ps = env.get_point(v, 's'); // ps = H(x, v) or H(u, v)

      // For every t_K, calculate f(x, v) 
      // tmp is every thing in f(x, v) except the prop involving t_K
      LatticePropagator tmp_D1Q1K(env.grid), tmp_D1Q2K(env.grid), tmp_D1Q1Kbar(env.grid), tmp_D1Q2Kbar(env.grid);  
      LatticePropagator tmp_D2Q1K(env.grid), tmp_D2Q2K(env.grid), tmp_D2Q1Kbar(env.grid), tmp_D2Q2Kbar(env.grid);  
      LatticePropagator tmp_sBar_d_D1(env.grid), tmp_sBar_d_D2(env.grid);  
      tmp_D1Q1K = Zero(); tmp_D1Q2K = Zero(); tmp_D1Q1Kbar = Zero(); tmp_D1Q2Kbar = Zero();
      tmp_D2Q1K = Zero(); tmp_D2Q2K = Zero(); tmp_D2Q1Kbar = Zero(); tmp_D2Q2Kbar = Zero();
      for(int mu=0; mu<4; ++mu) {
        tmp_D1Q1K += trace(gL[mu] * Lxx) *  (gL[mu] * pl); // tmp_D1Q1K = sum_mu tr(gL L(x,x)) gL L(x, v)
        tmp_D1Q2K += gL[mu] * Lxx * gL[mu] * pl; 
        tmp_D1Q1Kbar += trace(gL[mu] * Lxx) *  adj(pl) * gL[mu]; 
        tmp_D1Q2Kbar += adj(pl) * gL[mu] * Lxx * gL[mu]; 

        tmp_D2Q1K += trace(gL[mu] * Lxx) * adj(ps) * gL[mu]; 
        tmp_D2Q2K += adj(ps) * gL[mu] * Lxx *  gL[mu]; 
        tmp_D2Q1Kbar += trace(gL[mu] * Lxx) * (gL[mu] * ps); 
        tmp_D2Q2Kbar += gL[mu] * Lxx * gL[mu] * ps; 
      }
      tmp_sBar_d_D1 = g5 * pl; 
      tmp_sBar_d_D2 = adj(ps) * g5; 

      // map: t_K -> f(x, v) // For given v, the possible interval of tK is [tv - T/2 + tsep2 + 1, tv - tsep]
      map<int, LatticePropagatorSite> fx_D1Q1K, fx_D1Q2K, fx_D1Q1Kbar, fx_D1Q2Kbar; 
      map<int, LatticePropagatorSite> fx_D2Q1K, fx_D2Q2K, fx_D2Q1Kbar, fx_D2Q2Kbar; 
      map<int, LatticePropagatorSite> fx_sBar_d_D1, fx_sBar_d_D2; 
      for(int tK = v[3] - T/2 + tsep2 + 1; tK <= v[3] - tsep; ++tK) { 
        int original_tK = tK;
        if(tK<0) tK += T; // tK in [0, T]

        int lower_bound = tK + tsep3, upper_bound = tK + T/2 - 2; // upper_bound can be >= T
        Sum_Interval sum_interval(lower_bound, upper_bound, T);

        fx_D1Q1K[tK] = sum_interval(LatticePropagator( adj(ws[tK]) * tmp_D1Q1K )); // sum_x H(x, t_K)^\dagger sum_mu tr(gL L(x,x)) gL L(x, v)
        fx_D1Q2K[tK] = sum_interval(LatticePropagator( adj(ws[tK]) * tmp_D1Q2K )); 
        fx_D1Q1Kbar[tK] = sum_interval(LatticePropagator( tmp_D1Q1Kbar * ws[tK] )); 
        fx_D1Q2Kbar[tK] = sum_interval(LatticePropagator( tmp_D1Q2Kbar * ws[tK] )); 

        fx_D2Q1K[tK] = sum_interval(LatticePropagator( tmp_D2Q1K * wl[tK] )); 
        fx_D2Q2K[tK] = sum_interval(LatticePropagator( tmp_D2Q2K * wl[tK] )); 
        fx_D2Q1Kbar[tK] = sum_interval(LatticePropagator( adj(wl[tK]) * tmp_D2Q1Kbar)); 
        fx_D2Q2Kbar[tK] = sum_interval(LatticePropagator( adj(wl[tK]) * tmp_D2Q2Kbar)); 

        fx_sBar_d_D1[tK] = sum_interval(LatticePropagator( adj(ws[tK]) * tmp_sBar_d_D1));
        fx_sBar_d_D2[tK] = sum_interval(LatticePropagator( tmp_sBar_d_D2 * wl[tK]));

        tK = original_tK;
      }

      // combine f(x,v) with g(u, v)
      LatticePGG rst_D1Q1K(env.grid), rst_D1Q2K(env.grid), rst_D1Q1Kbar(env.grid), rst_D1Q2Kbar(env.grid); 
      LatticePGG rst_D2Q1K(env.grid), rst_D2Q2K(env.grid), rst_D2Q1Kbar(env.grid), rst_D2Q2Kbar(env.grid); 
      LatticePGG rst_sBar_d_D1(env.grid), rst_sBar_d_D2(env.grid); 
      vector<LatticePGG*> rst_vec = {&rst_D1Q1K, &rst_D1Q2K, &rst_D1Q1Kbar, &rst_D1Q2Kbar, 
                                     &rst_D2Q1K, &rst_D2Q2K, &rst_D2Q1Kbar, &rst_D2Q2Kbar, 
                                     &rst_sBar_d_D1, &rst_sBar_d_D2}; // for brevity of code
      for(auto rst: rst_vec) *rst = Zero();

      thread_for(ss, env.grid->lSites(), {
        Coordinate lcoor, gcoor;
        localIndexToLocalGlobalCoor(env.grid, ss, lcoor, gcoor);

        if(distance(v[3], gcoor[3], T) >= T/2 - tsep - tsep2) continue; // |t_u - t_v| must be less than T/2 - tsep - tsep2

        int tK = leftPoint(v[3], gcoor[3], T) - tsep; // t_wall = min(t_u, t_v) - t_sep
        if(tK < 0) tK += T; 
        int dist_v_wall = distance(v[3], tK, T); // distance from wall to v; always positive

        LatticePropagatorSite wall_to_u_L, wall_to_u_S, v_to_u_L, v_to_u_S; // L(u, t_K), L(u, v)
        peekLocalSite(wall_to_u_L, wl[tK], lcoor);    peekLocalSite(wall_to_u_S, ws[tK], lcoor);
        peekLocalSite(v_to_u_L, pl, lcoor);           peekLocalSite(v_to_u_S, ps, lcoor);

        LatticePGGSite rst_D1Q1K_site, rst_D1Q2K_site, rst_D1Q1Kbar_site, rst_D1Q2Kbar_site;
        LatticePGGSite rst_D2Q1K_site, rst_D2Q2K_site, rst_D2Q1Kbar_site, rst_D2Q2Kbar_site;
        LatticePGGSite rst_sBar_d_D1_site, rst_sBar_d_D2_site;
        vector<LatticePGGSite *> rst_site_vec = {&rst_D1Q1K_site, &rst_D1Q2K_site, &rst_D1Q1Kbar_site, &rst_D1Q2Kbar_site,
                                                 &rst_D2Q1K_site, &rst_D2Q2K_site, &rst_D2Q1Kbar_site, &rst_D2Q2Kbar_site, 
                                                 &rst_sBar_d_D1_site, &rst_sBar_d_D2_site}; // for brevity of code

        for(int mu=0; mu<4; ++mu) {
          for(int nu=0; nu<4; ++nu) {
            LatticePropagatorSite gu_D1 = gmu5[nu] * adj(v_to_u_L) * gmu5[mu] * wall_to_u_L; // g(u, v)

            rst_D1Q1K_site()()(mu, nu) = trace(fx_D1Q1K.at(tK) * gu_D1);
            rst_D1Q2K_site()()(mu, nu) = trace(fx_D1Q2K.at(tK) * gu_D1);
            rst_D1Q1Kbar_site()()(mu, nu) = trace(fx_D1Q1Kbar.at(tK) * adj(gu_D1));
            rst_D1Q2Kbar_site()()(mu, nu) = trace(fx_D1Q2Kbar.at(tK) * adj(gu_D1));
            rst_sBar_d_D1_site()()(mu, nu) = trace(fx_sBar_d_D1.at(tK) * gu_D1);

            LatticePropagatorSite gu_D2 = adj(wall_to_u_S) * gmu5[mu] * v_to_u_S * gmu5[nu]; // g(u, v)

            rst_D2Q1K_site()()(mu, nu) = trace(fx_D2Q1K.at(tK) * gu_D2);
            rst_D2Q2K_site()()(mu, nu) = trace(fx_D2Q2K.at(tK) * gu_D2);
            rst_D2Q1Kbar_site()()(mu, nu) = trace(fx_D2Q1Kbar.at(tK) * adj(gu_D2));
            rst_D2Q2Kbar_site()()(mu, nu) = trace(fx_D2Q2Kbar.at(tK) * adj(gu_D2));
            rst_sBar_d_D2_site()()(mu, nu) = trace(fx_sBar_d_D2.at(tK) * gu_D2);
          }
        }

        for(auto site: rst_site_vec) *site = (*site) * std::exp(env.M_K * dist_v_wall);
        for(int i=0; i<rst_vec.size(); ++i) pokeLocalSite(*rst_site_vec[i], *rst_vec[i], lcoor);
      });

      for(int i=0; i<rst_vec.size(); ++i) {
        for(int mu=0; mu<4; ++mu) *rst_vec[i] = Cshift(*rst_vec[i], mu, v[mu]); // shift v to origin
        *rst_vec_allsrc[i] += *rst_vec[i];
      }

    } // end of point source loop
    
    for(auto rst: rst_vec_allsrc) *rst = *rst * (1. / double(env.N_pt_src));

    writeScidac(rst_D1Q1K_allsrc, env.out_prefix + "/typeIII/D1Q1K." + to_string(traj));
    writeScidac(rst_D1Q2K_allsrc, env.out_prefix + "/typeIII/D1Q2K." + to_string(traj));
    writeScidac(rst_D1Q1Kbar_allsrc, env.out_prefix + "/typeIII/D1Q1Kbar." + to_string(traj));
    writeScidac(rst_D1Q2Kbar_allsrc, env.out_prefix + "/typeIII/D1Q2Kbar." + to_string(traj));

    writeScidac(rst_D2Q1K_allsrc, env.out_prefix + "/typeIII/D2Q1K." + to_string(traj));
    writeScidac(rst_D2Q2K_allsrc, env.out_prefix + "/typeIII/D2Q2K." + to_string(traj));
    writeScidac(rst_D2Q1Kbar_allsrc, env.out_prefix + "/typeIII/D2Q1Kbar." + to_string(traj));
    writeScidac(rst_D2Q2Kbar_allsrc, env.out_prefix + "/typeIII/D2Q2Kbar." + to_string(traj));

    writeScidac(rst_sBar_d_D1_allsrc, env.out_prefix + "/JJ_sBar_d_K/D1." + to_string(traj));
    writeScidac(rst_sBar_d_D2_allsrc, env.out_prefix + "/JJ_sBar_d_K/D2." + to_string(traj));

    // writeScidac(t1, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/KGG/KGG_typeIII." + to_string(traj));
  
  } // end of traj loop

  return 0;
}

