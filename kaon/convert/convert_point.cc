#include <headers/headers.h>
#include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

Coordinate gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  Coordinate mpi_coor = GridDefaultMpi();
  qlat::begin(&argc, &argv, qlat::Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor); 
  int T = grid->_fdimensions[Tdir];
  LatticePropagator prop(grid);

  // int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon
  int traj_start = 2400, traj_end = 2400, traj_sep = 100; // for 24ID, kaon
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  std::vector<std::vector<int>> xgs_l;
  std::vector<std::vector<int>> xgs_s;
  std::map<std::vector<int>, std::string> point_subdirs_l;
  std::map<std::vector<int>, std::string> point_subdirs_s;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    string point_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/luchang/point/prop-hvp ; results=" + std::to_string(traj)  + "/huge-data/prop-point-src";
    get_xgs(point_path, xgs_l, point_subdirs_l, 'l');
    get_xgs(point_path, xgs_s, point_subdirs_s, 's');

    // xgs_s.resize(1);
    for(auto x: xgs_s) {
      string path = point_subdirs_l[x]; // l quark wall source
      std::string out_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/point_l/" + std::to_string(traj) + "/" + coor2CSL(x);
      read_qlat_propagator_no_dist(prop, path); //  use no_dist
      writeScidac_prop_d2f(prop, out_path);

      path = point_subdirs_s[x]; // l quark wall source
      out_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/point_s/" + std::to_string(traj) + "/" + coor2CSL(x);
      read_qlat_propagator_no_dist(prop, path);
      writeScidac_prop_d2f(prop, out_path);
    }
  }

  qlat::end();

  return 0;
}
