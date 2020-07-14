
// On 8 nodes, Needs 16h for one trajectory (500 piont sources) // ~110s per point source

#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

accelerator_inline LatticePropagatorSite outerProduct(const LatticeFermionSite &f1, const LatticeFermionSite &f2) {
  LatticePropagatorSite rst;
  for(int s1=0; s1<4; s1++)
    for(int s2=0; s2<4; s2++)
      for(int c1=0; c1<3; c1++)
        for(int c2=0; c2<3; c2++)
            rst()(s1, s2)(c1, c2) = f1()(s1)(c1) * conjugate(f2()(s2)(c2));
            // rst()(s1, s2)(c1, c2) = f1()(s1)(c1) * std::conj(f2()(s2)(c2));
  return rst;
}

LatticePropagator outerProduct(const LatticeFermionSite &site, const LatticeFermion &lat) {

  LatticePropagator rst(lat.Grid());

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);

    LatticeFermionSite lat_s;
    peekLocalSite(lat_s, lat, lcoor);

    LatticePropagatorSite m = outerProduct(site, lat_s);

    pokeLocalSite(m, rst, lcoor);
  });

  return rst;
}



std::vector<int> gcoor({24, 24, 24, 64});


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  // int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  // FIXME: change those parameters
  int tsep = 16;

  Env env(gcoor, "24ID");
  // init_para(argc, argv, env);
  // env.N_pt_src = 1;  // FIXME: keep only one point
  env.N_pt_src = -1;  // Use all points

  const int T = env.grid->_fdimensions[3];

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    std::vector<LatticePropagator> wl = env.get_wall('l');
    std::vector<LatticePropagator> ws = env.get_wall('s');

    LatticeKGG rst_Q1_allsrc(env.grid), rst_Q2_allsrc(env.grid); 
    vector<LatticeKGG*> rst_vec_allsrc= {&rst_Q1_allsrc, &rst_Q2_allsrc}; 
    for(auto rst: rst_vec_allsrc) *rst = Zero();

    std::vector<LatticeFermionD> a2a_v = env.get_a2a('v');
    std::vector<LatticeFermionD> a2a_w = env.get_a2a('w');

    int num_pt_src = 0;
    if(env.N_pt_src != -1) env.xgs_s.resize(env.N_pt_src);
    for(const auto &v: env.xgs_s) {
      ++num_pt_src;

      LatticePropagator pl = env.get_point(v, 'l'); // pl = L(x, v) 

      int tK = v[3] - tsep;
      if(tK < 0) tK += T;

      vector<LatticePropagator> fx(4, env.grid);   // f(x) = gL * (L(x, tK) * H(x, tK)^dagger - H(x, tK) * L(x, tK)^dagger )
      LatticePropagator tmp = wl[tK] * adj(ws[tK]);
      tmp = tmp - adj(tmp);
      for(int rho=0; rho<4; ++rho) fx[rho] = gL[rho] * tmp;

      LatticePropagator tmp_Q1(env.grid), tmp_Q2(env.grid);
      tmp_Q1 = Zero(); tmp_Q2 = Zero();
      for(int rho=0; rho<4; ++rho) {
        tmp_Q1 += trace(fx[rho]) * adj(pl) * gL[rho];
        tmp_Q2 += adj(pl) * fx[rho] * gL[rho];
      }

      LatticePropagator huv_Q1(env.grid), huv_Q2(env.grid);
      huv_Q1 = Zero(); huv_Q2 = Zero();
      for(int i=0; i<a2a_v.size(); ++i) {

        if(i % 200 == 0) std::cout << GridLogMessage <<"Summing A2A modes " << i << std::endl;

        LatticeFermionSite gv_i_Q1, gv_i_Q2;
        gv_i_Q1 = sum(LatticeFermion(tmp_Q1 * a2a_v[i]));
        gv_i_Q2 = sum(LatticeFermion(tmp_Q2 * a2a_v[i]));

        huv_Q1 += outerProduct(gv_i_Q1, adj(a2a_w[i]));
        huv_Q2 += outerProduct(gv_i_Q2, adj(a2a_w[i]));

      }

      LatticeKGG rst_Q1(env.grid), rst_Q2(env.grid); 
      vector<LatticeKGG*> rst_vec = {&rst_Q1, &rst_Q2};
      for(auto rst: rst_vec) *rst = Zero();

      for(int mu=0; mu<4; ++mu) {
        for(int nu=0; nu<4; ++nu) {
          LatticeComplex tmp = trace(huv_Q1 * gmu[mu] * pl * gmu5[nu]);
          PokeIndex<LorentzIndex>(rst_Q1, tmp, mu, nu);

          tmp = trace(huv_Q2 * gmu[mu] * pl * gmu5[nu]);
          PokeIndex<LorentzIndex>(rst_Q2, tmp, mu, nu);
        }
      }
      // std::cout << "rst: " << std::endl;
      // print_grid_field_site(rst_Q1, {0,0,0,0});
      // print_grid_field_site(rst_Q2, {0,0,0,0});

      for(int i=0; i<rst_vec.size(); ++i) {
        *rst_vec[i] *= std::exp(env.M_K * tsep);   // add exp factor
        for(int mu=0; mu<4; ++mu) *rst_vec[i] = Cshift(*rst_vec[i], mu, v[mu]); // shift v to origin 
        *rst_vec_allsrc[i] += *rst_vec[i];
      }


    } // end of point source loop

    for(auto rst: rst_vec_allsrc) *rst = *rst * (1. / double(num_pt_src));

    writeScidac(rst_Q1_allsrc, env.out_prefix + "/typeII/Q1." + to_string(traj));
    writeScidac(rst_Q2_allsrc, env.out_prefix + "/typeII/Q2." + to_string(traj));

  } // end of traj loop

  std::cout << "Finished!" << std::endl;
  Grid_finalize();

  return 0;
}
