// #include <headers/headers.h>
#include "../kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


Coordinate gcoor({24, 24, 24, 64});


int main(int argc, char* argv[])
{
  zyd_init_Grid_Qlattice(argc, argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi()); 
  LatticePropagator prop(grid);

  int target_traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) { 
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, target_traj);
  }
  else {
    std::cout << "traj not specified; exiting" << std::endl;
    assert(0);
  }
  int traj_start = target_traj;
  int traj_end = target_traj;
  int traj_sep = 10; 

  // int traj_start = 2160, traj_end = 2160, traj_sep = 100; // for 24ID, kaon
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_num: " << traj_num << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  std::vector<std::vector<int>> xgs_l;
  std::map<std::vector<int>, std::string> point_subdirs_l;
  std::vector<std::vector<int>> xgs_s;
  std::map<std::vector<int>, std::string> point_subdirs_s;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    // string point_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/luchang/point/prop-hvp ; results=" + std::to_string(traj)  + "/huge-data/prop-point-src";
    string point_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_luchang_props/point/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
    get_xgs(point_path, xgs_l, point_subdirs_l, 'l');
    get_xgs(point_path, xgs_s, point_subdirs_s, 's');

    // xgs_s.resize(1);
    for(auto x: xgs_s) {
      string path = point_subdirs_l[x]; // l quark wall source
      std::string out_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/point_l/" + std::to_string(traj) + "/" + coor2str(x);
      // read_qlat_propagator_no_dist(prop, path); //  use no_dist
      read_qlat_propagator(prop, path);
      writeScidac_prop_d2f(prop, out_path);

      path = point_subdirs_s[x]; // l quark wall source
      out_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_my_props/point_s/" + std::to_string(traj) + "/" + coor2str(x);
      // read_qlat_propagator_no_dist(prop, path);
      read_qlat_propagator(prop, path);
      writeScidac_prop_d2f(prop, out_path);
    }
  }
  std::cout << "Number of point sources: " << xgs_s.size() << std::endl;
  std::cout << "Finished!" << std::endl;

  qlat::end();

  return 0;
}
