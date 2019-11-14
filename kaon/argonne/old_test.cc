#include <headers/headers.h>
#include "env.h"
#include "util.h"


using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});

template<typename vtype>
inline iColourMatrix<vtype> traceS(const iSpinColourMatrix<vtype> &p) {
  iColourMatrix<vtype> rst;
  rst()() = p()(0,0);
  for(int mu=1; mu<4; ++mu) rst()() += p()(mu,mu);
  return rst;
}

LatticeColourMatrix traceS(const LatticePropagator &p) {
  LatticeColourMatrix rst(p._grid);
  parallel_for(int ss=0; ss<p._grid->oSites(); ++ss)
    rst._odata[ss] = traceS(p._odata[ss]);
  return rst;
}

inline iMatrix<iScalar<iScalar<vComplex>>, 4> traceC(const iMatrix<iScalar<iMatrix<vComplex, 3> >, 4> &p) {
  iMatrix<iScalar<iScalar<vComplex>>, 4> rst;
  rst = zero;
  for(int mu=0; mu<4; ++mu)
    for(int nu=0; nu<4; ++nu)
      for(int c=0; c<3; ++c) {
        rst(mu, nu)()() += p(mu, nu)()(c, c);
      }

  return rst;
}


LatticeKGG traceC(const LatticeLorentzColour &p) {
  LatticeKGG rst(p._grid);
  parallel_for(int ss=0; ss<p._grid->oSites(); ++ss) {
    rst._odata[ss] = traceC(p._odata[ss]);
  }
  return rst;
}

using LatticeKGG = Lattice<iMatrix<iScalar<iScalar<vComplex> >, 4>>;
using LatticeLorentzColour = Lattice<iMatrix<iScalar<iMatrix<vComplex, 3> >, 4>>;

// TODO: 1. add exponential factor stemming from position of Kaon 2. For now, I only do it for one x  3. Do not need to calculate wall source on every site
LatticeKGG typeI(const Env &env, const std::vector<int> &v, int t_min) {

  LatticeKGG rst(env.grid);
  LatticeKGG rst_q1_avg(env.grid);
  LatticeKGG rst_q2_avg(env.grid);

    std::cout << GridLogMessage << std::endl;
    std::vector<LatticePropagator> wl = env.get_wall_l();
    std::vector<LatticePropagator> ws = env.get_wall_s();
    std::cout << GridLogMessage << std::endl;

  for(const auto &x: env.xgs_s) {
    std::cout << "point src x: " << x << std::endl;

    // the first loop J-Hw
    LatticePropagator pl = env.get_point_l(x); // pl = L(u, x)
    // LatticePGG loop1(env.grid);
    LatticeLorentzColour loop1(env.grid);
    for(int mu=0; mu<4; ++mu)
      for(int rho=0; rho<4; ++rho) {
        // LatticeComplex tmp(env.grid);
        // tmp = trace(gmu5[mu] * pl * gL[rho] * adj(pl));
        // pokeColour(loop1, tmp, mu, rho);
        LatticeColourMatrix tmp(env.grid);
        tmp = traceS(gmu5[mu] * pl * gL[rho] * adj(pl));
        pokeLorentz(loop1, tmp, mu, rho);
      }

    LatticeKGG loop1_t(env.grid);
    loop1_t = traceC(loop1);
    std::cout << loop1_t << std::endl;

    // the second loop J-Hw-K
    // LatticePGG loop2_1(env.grid), loop2_2(env.grid);
    LatticeLorentzColour loop2_1(env.grid), loop2_2(env.grid);
    LatticePropagator ps = env.get_point_s(x);

    int T = env.grid->_fdimensions[3];
    int t_wall = x[3] - t_min;
    if(t_wall < 0) t_wall += T; // if wall is on the left side of the current

    typename LatticePropagator::vector_object::scalar_object wall_to_x_l, wall_to_x_s;
    peekSite(wall_to_x_l, wl[t_wall], x);
    peekSite(wall_to_x_s, ws[t_wall], x);


    for(int rho=0; rho<4; ++rho)
      for(int nu=0; nu<4; ++nu) {
        // LatticeComplex tmp(env.grid);
        // tmp = trace(gL[rho] * adj(pl) * gmu5[nu] * wl[t_wall] * adj(wall_to_x_s)); // can be optimized by changing order and storing adj
        // pokeColour(loop2_1, tmp, rho, nu);
        // tmp = trace(gL[rho] * wall_to_x_l * adj(ws[t_wall]) * gmu5[nu] * ps);
        // pokeColour(loop2_2, tmp, rho, nu);
        LatticeColourMatrix tmp(env.grid);
        tmp = traceS(gL[rho] * adj(pl) * gmu5[nu] * wl[t_wall] * adj(wall_to_x_s)); // can be optimized by changing order and storing adj
        pokeLorentz(loop2_1, tmp, rho, nu);
        tmp = traceS(gL[rho] * wall_to_x_l * adj(ws[t_wall]) * gmu5[nu] * ps);
        pokeLorentz(loop2_2, tmp, rho, nu);
      }

    LatticeLorentzColour loop2(env.grid);
    loop2 = loop2_1 - loop2_2;

    // TBD: cross correlation -> 
    //

  // FFT theFFT((GridCartesian *)loop1._grid);
  // LatticeLorentzColour loop1_fft(loop1._grid), loop2_fft(loop1._grid);
  // theFFT.FFT_all_dim(loop1_fft, loop1, FFT::forward);
  // theFFT.FFT_all_dim(loop2_fft, loop2, FFT::forward);
  //
  // LatticeKGG rst_q1_fft(env.grid), rst_q2_fft(env.grid);
  // theFFT.FFT_all_dim(ret_u, ret_u_fft, FFT::backward);
  //
  // LatticeKGG rst_q1(env.grid), rst_q2(env.grid);
  //
  //   rst_q1 = traceC(loop1) * traceC(loop2);
  //   rst_q2 = traceC(loop1 * loop2);

    // LatticePGG loop2(env.grid);
    // loop2 = loop2_1 - loop2_2;
    // std::cout << loop2 << std::endl;

    // int t_wall = leftPoint(v[3], x[3], T) - t_min; // the first argument is the time to be shifted to 0
    // if(t_wall < 0) t_wall += T; // if wall is on the left side of the current
    // std::cout << "t wall" << t_wall << std::endl;
    //
    // LatticePropagator ps = env.get_point_s(x);
    // typename LatticePropagator::vector_object::scalar_object x_to_v_l, x_to_v_s, wall_to_v_l, wall_to_x_l, wall_to_v_s, wall_to_x_s;
    // peekSite(x_to_v_l, pl, v);
    // peekSite(x_to_v_s, ps, v);
    // peekSite(wall_to_v_l, wl[t_wall], v);
    // peekSite(wall_to_v_s, ws[t_wall], v);
    // peekSite(wall_to_x_l, wl[t_wall], x);
    // peekSite(wall_to_x_s, ws[t_wall], x);
    //
    // typename LatticePGG::vector_object::scalar_object l2_1, l2_2; // loop2 
    //
    // for(int rho=0; rho<4; ++rho)
    //   for(int nu=0; nu<4; ++nu) {
    //     l2_1()()(rho, nu) = trace(gL[rho] * adj(x_to_v_l) * gmu5[nu] * wall_to_v_l * adj(wall_to_x_s))()()();
    //     l2_2()()(rho, nu) = trace(gL[rho] * wall_to_x_l * adj(wall_to_v_s) * gmu5[nu] * x_to_v_s)()()();
    //   }

    // rst = (4. / 9.) * l1 * (l2_2 - l2_1);

    break;
  }

  std::cout << GridLogMessage << std::endl;
  return rst;
}

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  int traj_start = 2300, traj_end = 2300, traj_sep = 100; // for 24ID, kaon wall
  // int traj_start = 1000, traj_end = 1000, traj_sep = 100; // for 24ID, kaon point
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;


  int t_min = 16;
  Env env(gcoor, "24ID");
  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    env.setup_traj(traj);
    // cout << env.xgs_l << endl;;
    // cout << env.xgs_s << endl;;


    // /////////////////////

    std::vector<int> v {0,0,0,0};
    LatticeKGG t1 = typeI(env, v, t_min);
    // LatticePGG t1 = typeI(env, v, t_min);
    // std::cout << t1 << std::endl;

  }


  end();

  return 0;
}
