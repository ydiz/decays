// Did not do free field test for this diagram; instead, set EM factor to pion mass and compare with Pion form factor.
#include "../kaon.h"

using namespace std;
using namespace Grid;

// LatticeKGG calc_JJ_pion(const Coordinate &v, const std::vector<LatticePropagator> &wl, const LatticePropagator &pl, int tsep, double M_pi) {
//
//   GridBase *grid = pl.Grid();
//   const int T = grid->_fdimensions[3];
//
//   LatticeKGG rst(grid); 
//
//   std::vector<LatticePropagatorSite> wall_to_v(T);
//   for(int i=0; i<T; ++i)	peekSite(wall_to_v[i], wl[i], v); 
//
//   autoView(rst_v, rst, CpuWrite);
//   autoView(pl_v, pl, CpuRead);
//   using vobj = typename LatticePropagator::vector_object;
//   std::vector<LatticeView<vobj>> wl_v;
//   for(int t=0; t<T; ++t) {
//     autoView(view, wl[t], CpuRead)
//     wl_v.push_back(vierst_site(mu, nu)()()w);
//   }
//
//   // std::cout << GridLogMessage << "Before thread_for" << std::endl;
//   thread_for(ss, grid->lSites(), {    // iterate over point u
//     Coordinate lcoor, gcoor;
//     localIndexToLocalGlobalCoor(grid, ss, lcoor, gcoor);
//
//     int t_pi = leftPoint(v[3], gcoor[3], T) - tsep; // t_wall = min(t_u, t_v) - t_sep
//     if(t_pi < 0) t_pi += T; 
//     int dist_v_wall = distance(v[3], t_pi, T); // distance from wall to v; always positive
//
//     LatticePropagatorSite wall_to_u, v_to_u;
//     peekLocalSite(wall_to_u, wl_v[t_pi], lcoor);
//     peekLocalSite(v_to_u, pl_v, lcoor);
//
//     LatticeKGGSite rst_site;
//     for(int mu=0; mu<4; ++mu)
//       for(int nu=0rst_site(mu, nu)()(); nu<4; ++nu) {
//         rst_site(mu, nu)()() = 2. * real(trace(adj(wall_to_u) * gmu5[mu] * v_to_u * gmu[nu] * wall_to_v[t_pi]));
//       }
//
//     rst_site *= std::exp(M_pi * dist_v_wall);  // Note: here must be M_pi, not M_K !! 
//     pokeLocalSite(rst_site, rst_v, lcoor);
//
//   });
//   return rst;
// }

// return exp(m_pi (t_pi - t_v)) <J(u)J(v) pi(t_pi)> as a function of u; 
// This is later combined with half lattice Euv, so the region where t_u < t_v is not used
LatticeKGG calc_JJ_pion(const Coordinate &v, const std::vector<LatticePropagator> &wl, 
                        const LatticePropagator &pl, int tsep, double M_pi) 
{
  GridBase *grid = pl.Grid();
  const int T = grid->_fdimensions[3];

  int t_pi = (v[3] - tsep + T) % T; // t_wall = t_v - t_sep; we always assume that t_u >= t_v

  LatticePropagatorSite wall_to_v;
  peekSite(wall_to_v, wl[t_pi], v); 

  LatticeKGG rst(grid); rst = Zero();

  for(int mu=0; mu<4; ++mu) {
    for(int nu=0; nu<4; ++nu) {
      if(mu==3 || nu==3 || mu==nu) continue;  // Euv is 0 for these sites.
      LatticeComplex rst_mu_nu(grid);
      rst_mu_nu = 2. * real(trace(adj(wl[t_pi]) * gmu5[mu] * pl * gmu[nu] * wall_to_v)); // two diagrams that are complex conjugate to each other
      pokeLorentz(rst, rst_mu_nu, mu, nu);
    }
  }

  rst *= std::exp(M_pi * tsep);  // Note: here must be M_pi, not M_K !! 

  return rst;
}

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  // zyd_init_Grid_Qlattice(argc, argv);

#ifdef CUTH_FREE_FIELD
  vector<int> t_seps = {1,2,3}; 
  int max_uv_sep = 2;
  Env env("FreeField_8nt8");
#else
  vector<int> t_seps = {8, 10, 12, 14, 16};  
  // int max_uv_sep = 16;
  int max_uv_sep = 32;
  Env env("24ID");
#endif

  const int T = env.grid->_fdimensions[3];

  int N_t_seps = t_seps.size();
  std::cout << "t_seps: " << t_seps << std::endl;

  int target_traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, target_traj);
  }
  else {
    std::cout << "--traj not specified; exiting" << std::endl;
    exit(0);
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

  
  // env.N_pt_src = 1;  // FIXME: keep only one point
  env.N_pt_src = -1;  // Use all points

  int traj_idx = 0;
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    vector<vector<Complex>> table2d(N_t_seps); // table_rst[tsep][ut]; ut=0 is always the position of v
    for(auto &x: table2d) x.resize(T, 0.);

    std::vector<LatticePropagator> wl = env.get_wall('l');

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_l.resize(env.N_pt_src);

    Complex rst_allsrc = 0.;
    for(const auto &v: env.xgs_l) {
      ++num_pt_src;
      std::cout << "num_pt_src: " << num_pt_src << std::endl;

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(u, v)

      LatticeKGG Euv(env.grid);
      EM_factor_half_lattice(Euv, v, env.M_K, max_uv_sep);
      // EM_factor(Euv, v, env.M_K, max_uv_sep); // Must use kaon mass
      // EM_factor(Euv, v, env.M_pi, max_uv_sep); // FIXME: for test; 

      for(int t_sep_idx=0; t_sep_idx<t_seps.size(); ++t_sep_idx) { // iterate through all possible separation between eta and kaon
        int t_sep = t_seps[t_sep_idx];

        LatticeKGG JJ_pi = calc_JJ_pion(v, wl, pl, t_sep, env.M_pi);

        LatticeComplex rst(env.grid); rst = Zero();
        for(int mu=0; mu<4; ++mu) { 
          for(int nu=0; nu<4; ++nu) { 
            if(mu==nu || mu==3 || nu == 3) continue; // Under these conditions, E_{mu,nu} = 0
            LatticeComplex E_munu = peekLorentz(Euv, mu, nu);
            LatticeComplex H_munu = peekLorentz(JJ_pi, mu, nu);
            rst += E_munu * H_munu;
          }
        }

        vector<LatticeComplexSite> rst_ut;
        sliceSum(rst, rst_ut, Tdir);

        for(int ut=0; ut<T; ++ut) table2d[t_sep_idx][ut] += rst_ut[(ut+v[3])%T]()()();
      }

    } // end of point source loop
    
    std::cout << GridLogMessage << "Number of point sources: " << num_pt_src << std::endl;
    for(int t_sep_idx=0; t_sep_idx<t_seps.size(); ++t_sep_idx)  {
      for(int ut=0; ut<T; ++ut) {
        table2d[t_sep_idx][ut] /= double(num_pt_src);
      }
    }
    std::cout << "traj [" << traj << "] JJ_pion: " << table2d << std::endl;

    ++traj_idx;
  } // end of traj loop


  // // For Test: Change EM factor to pion mass!! ; then, the following amplitude can be compared with ABJ prediction
  // double hadron_coef = 1./ (3 * std::sqrt(2)) * env.Z_V * env.Z_V * 2. * env.M_pi / env.N_pi;
  // vector<vector<Complex>> pion_form_factors = table_rst;
  // for(auto &x: pion_form_factors) {
  //   for(auto &y: x) {
  //     y *= hadron_coef;
  //   }
  // }
  // std::cout << "Form Factor with different tsep: " << pion_form_factors << std::endl;
  // std::cout << "(If use pion mass) ABJ prediction is 0.2724 GeV^-1" << std::endl;

  std::cout << GridLogMessage << "Finished!" << std::endl;
  Grid_finalize();
  return 0;
}

