#include "kaon.h"

using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

std::vector<int> gcoor({24, 24, 24, 64});

void test_type2_q1(Env &env) {

  std::cout << "number of point sources: " << env.xgs_s.size() << std::endl;
  env.xgs_s.resize(env.num_points);
  std::cout << "number of point sources using: " << env.xgs_s.size() << std::endl;

  const int T = env.grid->_fdimensions[3];

  std::mt19937 eng(0);
  std::uniform_int_distribution<int> dist(0, env.xgs_s.size()-1);
  auto rand = [&eng, &dist]() {return dist(eng);};

  for(int i = 0; i<env.num_pairs; ++i) {

    std::vector<int> u = env.xgs_s[rand()];
    std::vector<int> v = env.xgs_s[rand()];
    while(abs_distance(u[3], v[3], T) > T/4) v = env.xgs_s[rand()]; // make sure |t_u - t_v| <= T/4

    std::vector<int> r(4);
    for(int i=0; i<4; ++i) r[i] = v[i] - u[i];
    std::cout << "point src u: " << u << " point src v: " << v << " r: " << r << std::endl;

    int t_wall = left_time(u[3], v[3], T) - env.T_wall_typeII;
    if(t_wall < 0) t_wall += T;
    std::cout << "t_wall: " << t_wall << std::endl;

    LatticePropagator pl_u = env.get_point(u, 'l'); // light quark propagator with source u
    LatticePropagator pl_v = env.get_point(v, 'l'); // light quark propagator with source v
    LatticePropagator wl = env.get_wall(t_wall, 'l');
    LatticePropagator ws = env.get_wall(t_wall, 's');
    // conjugateU(pl_u); conjugateU(pl_v); conjugateU(wl); conjugateU(ws);
    // std::cout << "have done complex conjugate to U" << std::endl;

    std::vector<LatticeComplex> loop2(4, env.grid);
    LatticePropagator tmp = wl * adj(ws);
    for(int rho=0; rho<4; ++rho) {
      loop2[rho] = trace(gL[rho] * tmp);
    }
    std::cout << "loop2: {1,1,1,1}" << std::endl;
    print_grid_field_site(loop2[0], {1,1,1,1});
    print_grid_field_site(loop2[1], {1,1,1,1});

    typename LatticePropagator::vector_object::scalar_object Luv; // L(u, v)
    peekSite(Luv, pl_v, u);

    LatticeKGG rst_q1_uvx(env.grid); // for fixed u and v, calculate for every x
    for(int mu=0; mu<4; ++mu) {
      for(int nu=0; nu<4; ++nu) {
        LatticeComplex rst_mu_nu(env.grid);
        rst_mu_nu = Zero();

        LatticePropagator tmp(env.grid);
        tmp = pl_u * (gmu[mu] * Luv * gmu5[nu]) * adj(pl_v);
        for(int rho=0; rho<4; ++rho) {
          rst_mu_nu = rst_mu_nu + trace(gL[rho] * tmp) * loop2[rho];
        }
        pokeLorentz(rst_q1_uvx, rst_mu_nu, mu, nu);
      }
    }
    std::cout << "rst_q1_uvx: {1,1,1,1}" << std::endl;
    print_grid_field_site(rst_q1_uvx, {1,1,1,1});

    // If Hw is too close to the wall/on the wrong side of wall, set to 0
    int right_min = 5;
    LatticeComplex mask(env.grid);
    zero_mask(mask, t_wall, right_min);
    rst_q1_uvx = mask * rst_q1_uvx;

    typename LatticeKGG::vector_object::scalar_object rst_q1_uv; // L(u, v)
    rst_q1_uv = sum(rst_q1_uvx); // sum over x

    int left_dist = left_distance(u[3], t_wall, T);
    std::cout << "left_dist(distance between u and wall): " << left_dist << std::endl;
    rst_q1_uv = rst_q1_uv * std::exp(env.M_K * left_dist); // translation factor

    rst_q1_uv = rst_q1_uv * (8. / 9.); // Do not forget the coefficient

    std::cout << "rst_q1_uv: " << std::endl;
    std::cout << rst_q1_uv << std::endl;

    assert(0);
  }

}


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  // int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon wall
  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  Env env(gcoor, "24ID");
  init_para(argc, argv, env);

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    test_type2_q1(env);
    // LatticeKGG t1 = test_type2_q1(env);
    // // LatticeKGG t1 = typeI(env, t_min);
    // writeScidac(t1, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/KGG/KGG_typeII." + to_string(traj));

    // print_grid_field_site(t1, {1,1,1,1});
    // LatticePGG t1 = typeI(env, v, t_min);
  }

  return 0;
}
