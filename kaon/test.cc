#include <headers/headers.h>
#include "env.h"

using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});

template<typename vtype>
iColourMatrix traceS(const iSpinColourMatrix<vtype> &p) {
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

// TODO: 1. add exponential factor stemming from position of Kaon 2. For now, I only do it for one x  3. Do not need to calculate wall source on every site
LatticePGG typeI(const Env &env, const std::vector<int> &v, int t_min) {

  LatticePGG rst(env.grid);

    std::cout << GridLogMessage << std::endl;
    std::vector<LatticePropagator> wl = env.get_wall_l();
    std::vector<LatticePropagator> ws = env.get_wall_s();
    std::cout << GridLogMessage << std::endl;

  for(const auto &x: env.xgs_l) {
    std::cout << "point src x: " << x << std::endl;

    // the first loop J-Hw
    LatticePropagator pl = env.get_point_l(x); // pl = L(u, x)
    LatticePGG loop1(env.grid);
    for(int mu=0; mu<4; ++mu)
      for(int rho=0; rho<4; ++rho) {
        LatticeComplex tmp(env.grid);
        tmp = trace(gmu5[mu] * pl * gL[rho] * adj(pl));
        pokeColour(loop1, tmp, mu, rho);
      }

    // the second loop J-Hw-K
    LatticePGG loop2(env.grid);
    int T = env.grid->_fdimensions[3];
    int t_wall = x[3] - T;

    // parallel_for(int ss=0; ss<ret._grid->lSites(); ss++){
    //
    //   std::vector<int> lcoor, gcoor;
    //   localIndexToLocalGlobalCoor(ret._grid, ss, lcoor, gcoor);
    //   std::vector<int> &v = gcoor;
    //   int t_wall = leftPoint(v[3], x[3], T) - t_min;
    // if(t_wall < 0) t_wall += T; // if wall is on the left side of the current
    //
    // typename LatticePropagator::vector_object::scalar_object x_to_v_l, x_to_v_s, wall_to_v_l, wall_to_x_l, wall_to_v_s, wall_to_x_s;
    // // peekLocalSite(wall_to_xp, wall_props[t_wall], lcoor);
    // peekLocalSite(x_to_v_l, pl, lcoor);
    // peekLocalSite(x_to_v_s, ps, v);
    // peekLocalSite(wall_to_v_l, wl[t_wall], v);
    // peekLocalSite(wall_to_v_s, ws[t_wall], v);
    // peekLocalSite(wall_to_x_l, wl[t_wall], x);
    // peekLocalSite(wall_to_x_s, ws[t_wall], x);
    // }

    int t_wall = leftPoint(v[3], x[3], T) - t_min; // the first argument is the time to be shifted to 0
    if(t_wall < 0) t_wall += T; // if wall is on the left side of the current
    std::cout << "t wall" << t_wall << std::endl;

    LatticePropagator ps = env.get_point_s(x);
    typename LatticePropagator::vector_object::scalar_object x_to_v_l, x_to_v_s, wall_to_v_l, wall_to_x_l, wall_to_v_s, wall_to_x_s;
    peekSite(x_to_v_l, pl, v);
    peekSite(x_to_v_s, ps, v);
    peekSite(wall_to_v_l, wl[t_wall], v);
    peekSite(wall_to_v_s, ws[t_wall], v);
    peekSite(wall_to_x_l, wl[t_wall], x);
    peekSite(wall_to_x_s, ws[t_wall], x);

    typename LatticePGG::vector_object::scalar_object l2_1, l2_2; // loop2 

    for(int rho=0; rho<4; ++rho)
      for(int nu=0; nu<4; ++nu) {
        l2_1()()(rho, nu) = trace(gL[rho] * adj(x_to_v_l) * gmu5[nu] * wall_to_v_l * adj(wall_to_x_s))()()();
        l2_2()()(rho, nu) = trace(gL[rho] * wall_to_x_l * adj(wall_to_v_s) * gmu5[nu] * x_to_v_s)()()();
      }

    rst = (4. / 9.) * l1 * (l2_2 - l2_1);

    break;
  }

  for(int mu=0; mu<4; ++mu) {
    if(v[mu]!=0) rst = Cshift(rst, mu, v[mu]); // shift v to 0
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
    LatticePGG t1 = typeI(env, v, t_min);
    std::cout << t1 << std::endl;

  }


  end();

  return 0;
}
