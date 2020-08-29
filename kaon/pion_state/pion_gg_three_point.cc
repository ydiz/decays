#include "../kaon.h"
#include "../../amplitude/form_factor.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// "hadronic" should be e^{M_pi |t_pi|} <J(u) J(0) pi(t_pi) >, rather than <J(-w/2) J(w/2) |pi>
void measure_pion_form_factor(const LatticeKGG &hadronic, const Env &env) {
      static LatticeKGG lep(env.grid);
      static bool initialized = false;
      if(!initialized) {
        form_factor_integrand(lep, env.M_pi);
        initialized = true;
      }

      double hadron_coef, lep_coef;
      lep_coef = 2. / std::pow(env.M_pi, 4);
      hadron_coef = 1./ (3 * std::sqrt(2)) * env.Z_V * env.Z_V * 2. * env.M_pi / env.N_pi;
      std::string cutoff_type = "time";

      std::vector<double> amplitude = form_factor(hadronic, lep, hadron_coef, lep_coef, cutoff_type);

      const int T = env.grid->_fdimensions[3];
      amplitude.resize(T/4 + 1);   // keep [0, T/4]
      std::cout << "Form Factor vs cutoff: " << amplitude << std::endl;
      std::cout << "ABJ prediction is 0.2724 GeV^-1" << std::endl;
}


// the leptonic part must be M_K
// the lep_coef must be M_K
// However, the hadron_coef must be M_pi and N_pi
double combine_kaon_form_factor(const LatticeKGG &hadronic, const Env &env) {
      static LatticeKGG lep(env.grid);
      static bool initialized = false;
      if(!initialized) {
        form_factor_integrand(lep, env.M_K);
        initialized = true;
      }

      double hadron_coef, lep_coef;
      lep_coef = 2. / std::pow(env.M_K, 4);
      hadron_coef = env.Z_V * env.Z_V * 2. * env.M_pi / env.N_pi;
      std::string cutoff_type = "time";

      std::vector<double> amplitude = form_factor(hadronic, lep, hadron_coef, lep_coef, cutoff_type);

      const int T = env.grid->_fdimensions[3];
      amplitude.resize(T/4 + 1);   // keep [0, T/4]
      std::cout << "Luv(u, v=0, M_K) * <J(u)J(0)|pi> vs cutoff: " << amplitude << std::endl;
      return amplitude.back();
}








int main(int argc, char* argv[])
{
  // Grid_init(&argc, &argv);
  zyd_init_Grid_Qlattice(argc, argv);

  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  int tsep = 10;
  std::cout << "t_sep: " << tsep << std::endl;

  Env env("24ID");
  // init_para(argc, argv, env);
  // env.N_pt_src = 1;  // FIXME: keep only one point
  env.N_pt_src = -1;  // Use all points

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    // std::cout << "wl" << std::endl;
    // print_grid_field_site(wl[0], {1,2,3,4});

    // save average of all point sources
    LatticeKGG rst_allsrc(env.grid); 
    rst_allsrc = Zero();

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);
    for(const auto &v: env.xgs_l) {
      ++num_pt_src;
      std::cout << "num_pt_src: " << num_pt_src << std::endl;

      std::vector<LatticePropagatorSite> wall_to_v(T);
      for(int i=0; i<T; ++i)	peekSite(wall_to_v[i], wl[i], v); 

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(u, v)

      LatticeKGG rst(env.grid); 

      // std::cout << GridLogMessage << "Before thread_for" << std::endl;
      thread_for(ss, env.grid->lSites(), {    // iterate over point u
        Coordinate lcoor, gcoor;
        localIndexToLocalGlobalCoor(env.grid, ss, lcoor, gcoor);

        int t_pi = leftPoint(v[3], gcoor[3], T) - tsep; // t_wall = min(t_u, t_v) - t_sep
        if(t_pi < 0) t_pi += T; 
        int dist_v_wall = distance(v[3], t_pi, T); // distance from wall to v; always positive

        LatticePropagatorSite wall_to_u, v_to_u;
        peekLocalSite(wall_to_u, wl[t_pi], lcoor);
        peekLocalSite(v_to_u, pl, lcoor);

        LatticeKGGSite rst_site;
        for(int mu=0; mu<4; ++mu)
          for(int nu=0; nu<4; ++nu) {
            // rst_site(mu, nu)()() = trace(adj(wall_to_u) * gmu5[mu] * v_to_u * gmu[nu] * wall_to_v[t_pi] );
            rst_site(mu, nu)()() = 2. * real(trace(adj(wall_to_u) * gmu5[mu] * v_to_u * gmu[nu] * wall_to_v[t_pi]));
          }

        rst_site *= std::exp(env.M_pi * dist_v_wall);
        pokeLocalSite(rst_site, rst, lcoor);

      });

      for(int mu=0; mu<4; ++mu) rst = Cshift(rst, mu, v[mu]); // shift v to origin
      rst_allsrc += rst;

      // measure_pion_form_factor(rst, env);
      combine_kaon_form_factor(rst, env);

    } // end of point source loop
    
    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;
    rst_allsrc *= 1. / double(num_pt_src);

    writeScidac(rst_allsrc, env.out_prefix + "/pion_gg/" + to_string(traj));

    // measure_pion_form_factor(rst_allsrc, env);
    combine_kaon_form_factor(rst_allsrc, env);
  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}

