#include <headers/headers.h>
#include "env.h"
#include "util.h"
#include "typeI.h"



using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});



LatticeKGG test_t1_q1(Env &env, int t_min) {

  LatticeKGG rst_q1_avg(env.grid), rst_q2_avg(env.grid);
  rst_q1_avg = zero; rst_q2_avg = zero;

  env.xgs_s.erase(env.xgs_s.begin()+1, env.xgs_s.end()); // FIXME: keep only one point
  for(const auto &x: env.xgs_s) {
    std::cout << "point src x: " << x << std::endl;

    // the first loop J-Hw
    LatticePropagator pl = env.get_point_l(x); // pl = L(u, x)

    LatticeKGG loop1(env.grid);
    for(int mu=0; mu<4; ++mu)
      for(int rho=0; rho<4; ++rho) {
        LatticeComplex tmp(env.grid);
        tmp = trace(gL[rho] * adj(pl) * gmu5[mu] * pl);
        // print_grid_field_site(tmp, {1,2,3,4,});
        pokeLorentz(loop1, tmp, mu, rho);
      }
    print_grid_field_site(loop1, {1,2,3,4});

    // the second loop J-Hw-K
    LatticeKGG loop2(env.grid);
    LatticePropagator ps = env.get_point_s(x);

    int T = env.grid->_fdimensions[3];
    int t_wall = x[3] - t_min;
    if(t_wall < 0) t_wall += T;

    LatticePropagator wl = env.get_wall_l(t_wall);
    LatticePropagator ws = env.get_wall_s(t_wall);

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
    print_grid_field_site(loop2, {1,2,3,4});

    // // before doing cross correlation, set sites on the left of the wall and close to the wall to be zero
    // LatticeComplex mask(env.grid);
    // zero_mask(mask, t_wall);
    // loop1 = loop1 * mask;
    // loop2 = loop2 * mask;


    std::cout << GridLogMessage << "before fft" << std::endl;
    FFT theFFT((GridCartesian *)loop1._grid);
    // LatticeLorentzColour loop1_fft(loop1._grid), loop2_fft(loop1._grid);
    LatticeKGG loop1_fft(loop1._grid), loop2_fft(loop1._grid);
    theFFT.FFT_all_dim(loop1_fft, loop1, FFT::forward);
    theFFT.FFT_all_dim(loop2_fft, loop2, FFT::forward);

    LatticeKGG rst_q1_fft(env.grid);
    rst_q1_fft = conjugate(loop1_fft) * loop2_fft;

    LatticeKGG rst_q1(env.grid);
    theFFT.FFT_all_dim(rst_q1, rst_q1_fft, FFT::backward);
    // do not need to divide it by V; cancel by the V factor coming from summation over x (position of Hw)

    std::cout << GridLogMessage << "after fft" << std::endl;

    print_grid_field_site(rst_q1, {1,2,3,4});

    rst_q1 = rst_q1 * std::exp(- env.M_K * t_wall);
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
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  // int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon wall
  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  // int traj_start = 1000, traj_end = 1000, traj_sep = 100; // for 24ID, kaon point
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;


  int t_min = 20;
  Env env(gcoor, "24ID");
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);

    // LatticeKGG t1 = typeI(env, t_min);
    // LatticeKGG t1 = test_t1_q1(env, t_min);
    // writeScidac(t1, "/projects/CSC249ADSE03/yidizhao/KGG_config/24ID/typeI/KGG_typeI." + to_string(traj));
    // writeScidac(t1, "./KGG_typeI." + to_string(traj));
    // writeScidac(t1, "./KGG_typeI." + to_string(traj) + ".no-mask");
    LatticeKGG t1 = typeI(env, t_min);
    writeScidac(t1, "./KGG_typeI." + to_string(traj) + ".typeI_fnc");

    // std::cout << t1 << std::endl;
    // print_grid_field_site(t1, {1,1,1,1});
    // print_grid_field_site(t1, {1,2,3,4});
    // LatticePGG t1 = typeI(env, v, t_min);
    // std::cout << t1 << std::endl;

  }


  end();

  return 0;
}
