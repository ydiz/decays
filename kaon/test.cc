#include "kaon.h"

using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});

LatticeKGG test_t1_q1(Env &env, int t_min, int T_u) {

  LatticeKGG rst_q1_avg(env.grid), rst_q2_avg(env.grid);
  rst_q1_avg = Zero(); rst_q2_avg = Zero();

  // env.xgs_s.erase(env.xgs_s.begin()+1, env.xgs_s.end()); // FIXME: keep only one point
  std::cout << "number of point sources: " << env.xgs_s.size() << std::endl;
  // env.xgs_s.resize(100);
  env.xgs_s.resize(30);
  for(const auto &x: env.xgs_s) {
    std::cout << "point src x: " << x << std::endl;

    // the first loop J-Hw
    LatticeKGG loop1(env.grid);
    LatticePropagator pl = env.get_point(x, 'l'); // pl = L(u, x)
    for(int mu=0; mu<4; ++mu)
      for(int rho=0; rho<4; ++rho) {
        LatticeComplex tmp(env.grid);
        tmp = trace(gL[rho] * adj(pl) * gmu5[mu] * pl);
        // print_grid_field_site(tmp, {1,2,3,4,});
        pokeLorentz(loop1, tmp, mu, rho);
      }
    // print_grid_field_site(loop1, {1,2,3,4});

    // the second loop J-Hw-K
    LatticeKGG loop2(env.grid);
    LatticePropagator ps = env.get_point(x, 's');

    int T = env.grid->_fdimensions[3];
    int t_wall = x[3] - t_min;
    if(t_wall < 0) t_wall += T;
    std::cout << "t_wall: " << t_wall << std::endl;

    LatticeComplex exp_factor(env.grid);
    expUMinusTwall(exp_factor, t_wall, x[3], T_u, env.M_K);
    // std::cout << exp_factor << std::endl;
    // assert(0);
    loop1 = loop1 * exp_factor;

    LatticePropagator wl = env.get_wall(t_wall, 'l');
    LatticePropagator ws = env.get_wall(t_wall, 's');

    typename LatticePropagator::vector_object::scalar_object wall_to_x_s;
    // peekSite(wall_to_x_s, ws[t_wall], x);
    peekSite(wall_to_x_s, ws, x);

    for(int rho=0; rho<4; ++rho)
      for(int nu=0; nu<4; ++nu) {
        LatticeComplex tmp(env.grid);
        // tmp = trace(gL[rho] * adj(pl) * gmu5[nu] * wl[t_wall] * adj(wall_to_x_s));
        tmp = trace(gL[rho] * adj(pl) * gmu5[nu] * wl * adj(wall_to_x_s));
        pokeLorentz(loop2, tmp, rho, nu);
      }
    // print_grid_field_site(loop2, {1,2,3,4});

    // // before doing cross correlation, set sites on the left of the wall and close to the wall to be zero
    // LatticeComplex mask(env.grid);
    // zero_mask(mask, t_wall);
    // loop1 = loop1 * mask;
    // loop2 = loop2 * mask;


    std::cout << GridLogMessage << "before fft" << std::endl;
    FFT theFFT((GridCartesian *)loop1.Grid());
    // LatticeLorentzColour loop1_fft(loop1._grid), loop2_fft(loop1._grid);
    LatticeKGG loop1_fft(loop1.Grid()), loop2_fft(loop1.Grid());
    theFFT.FFT_all_dim(loop1_fft, loop1, FFT::forward);
    theFFT.FFT_all_dim(loop2_fft, loop2, FFT::forward);

    LatticeKGG rst_q1_fft(env.grid);
    rst_q1_fft = conjugate(loop1_fft) * loop2_fft;

    LatticeKGG rst_q1(env.grid);
    theFFT.FFT_all_dim(rst_q1, rst_q1_fft, FFT::backward);
    // do not need to divide it by V; cancel by the V factor coming from summation over x (position of Hw)
    std::cout << GridLogMessage << "after fft" << std::endl;

    // print_grid_field_site(rst_q1, {1,2,3,4});

    // rst_q1 = rst_q1 * std::exp(- env.M_K * t_wall); // FIXME:
    print_grid_field_site(rst_q1, {1,2,3,4});
    rst_q1_avg += rst_q1;
  }


  LatticeKGG rst(env.grid);
  rst = env.wilson_c1 * rst_q1_avg * (1. / (double) env.xgs_s.size());
  rst = rst * (- 4. / 9.); // Do not forget the coefficient
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
  int T_wall = 20;
  int T_u = 4;
  Env env(gcoor, "24ID");
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    LatticeKGG t1 = test_t1_q1(env, T_wall, T_u);
    // LatticeKGG t1 = typeI(env, t_min);
    writeScidac(t1, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/KGG/KGG_typeI." + to_string(traj));

    // print_grid_field_site(t1, {1,1,1,1});
    // LatticePGG t1 = typeI(env, v, t_min);
  }

  return 0;
}
