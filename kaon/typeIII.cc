#include "kaon.h"

using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

std::vector<int> gcoor({24, 24, 24, 64});

LatticeKGG test_type3_q1(Env &env) {

  LatticeKGG rst_q1_avg(env.grid);
  rst_q1_avg = Zero();

  std::cout << "number of point sources: " << env.xgs_s.size() << std::endl;
  env.xgs_s.resize(env.num_points);
  std::cout << "number of point sources using: " << env.xgs_s.size() << std::endl;

  const int T = env.grid->_fdimensions[3];

  for(const auto &x: env.xgs_s) {
    std::cout << "point src x: " << x << std::endl;

    LatticeKGG rst_q1(env.grid);

    int t_wall = x[3] - env.T_wall;
    if(t_wall < 0) t_wall += T;
    std::cout << "t_wall: " << t_wall << std::endl;

    LatticePropagator ps = env.get_point(x, 's'); // strange quark propagator
    LatticePropagator pl = env.get_point(x, 'l'); // light quark propagator
    LatticePropagator wl = env.get_wall(t_wall, 'l');
    LatticePropagator ws = env.get_wall(t_wall, 's');

    LatticeComplex exp_factor(env.grid);
    expUMinusTwall(exp_factor, t_wall, x[3], env.T_u, env.M_K);

    std::vector<LatticePropagator> Fux(4, env.grid); // F_mu(u, x)
    typename LatticePropagator::vector_object::scalar_object Lxx; // L(x, x)
    peekSite(Lxx, ws, x);
    for(int mu=0; mu<4; ++mu) {
      Fux[mu] = Zero();
      LatticePropagator tmp = adj(pl) * gmu5[mu] * wl;
      for(int rho=0; rho<4; ++rho) {
        Fux[mu] += trace(gL[rho] * Lxx) * (gL[rho] * tmp);
      }
      Fux[mu] = Fux[mu] * exp_factor;
    }
    // print_grid_field_site(loop1, {1,2,3,4});

    // the second loop J-Hw-K
    std::vector<LatticePropagator> Gvx(4, env.grid); // G_nu(v, x)
    for(int nu=0; nu<4; ++nu) {
      Gvx[nu] = adj(ws) * gmu5[nu] * ps;
    }

    std::cout << GridLogMessage << "before fft" << std::endl;
    FFT theFFT((GridCartesian *)env.grid);
    std::vector<LatticePropagator> Fux_fft(4, env.grid), Gvx_fft(4, env.grid); // F_mu(u, x)
    for(int mu=0; mu<4; ++mu) {
      theFFT.FFT_all_dim(Fux_fft[mu], Fux[mu], FFT::forward);
      theFFT.FFT_all_dim(Gvx_fft[mu], Gvx[mu], FFT::forward);
    }
    for(int mu=0; mu<4; ++mu) {
      for(int nu=0; nu<4; ++nu) {
        LatticePropagator tmp_fft = conjugate(Fux_fft[mu]) * Gvx_fft[nu];
        LatticePropagator tmp(env.grid);
        theFFT.FFT_all_dim(tmp, tmp_fft, FFT::backward);

        LatticeComplex rst_mu_nu = trace(tmp);
        pokeLorentz(rst_q1, rst_mu_nu, mu, nu);
      }
    }

    std::cout << GridLogMessage << "after fft" << std::endl;

    // do not need to divide it by V; cancel by the V factor coming from summation over x (position of Hw)
    rst_q1 = rst_q1 * (T / double(2 * env.T_u + 1)); //FIXME:  // number of pairs to sum over is not V T. It is V * (2 * T_u + 1) 

    print_grid_field_site(rst_q1, {1,2,3,4});
    rst_q1_avg += rst_q1;
  }


  LatticeKGG rst(env.grid);
  rst = env.wilson_c1 * rst_q1_avg * (1. / (double) env.xgs_s.size());
  rst = rst * (- 2. / 9.); // Do not forget the coefficient
  return rst;
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

  // int t_min = 20;
  // int t_min = 12;
  // int T_wall = 20;
  // int T_u = 4;
  Env env(gcoor, "24ID");
  init_para(argc, argv, env);

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    LatticeKGG t1 = test_type3_q1(env);
    // LatticeKGG t1 = typeI(env, t_min);
    writeScidac(t1, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/KGG/KGG_typeIII." + to_string(traj));

    // print_grid_field_site(t1, {1,1,1,1});
    // LatticePGG t1 = typeI(env, v, t_min);
  }

  return 0;
}
