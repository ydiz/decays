#include <headers/headers.h>
// #include "kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

Coordinate gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::cout << "after init" << std::endl;
  Coordinate mpi_coor = GridDefaultMpi();
  qlat::begin(&argc, &argv, qlat::Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplexD::Nsimd()), mpi_coor); 
  // GridCartesian * grid_f = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplexF::Nsimd()), mpi_coor); 
  int T = grid->_fdimensions[Tdir];
  LatticePropagator prop(grid);
  // LatticePropagatorF prop_f(grid_f);

  int traj_start = 2300, traj_end = 2400, traj_sep = 100; // for 24ID, kaon
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {

    std::string gt_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/luchang/gauge_transform_and_wall_l/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
    // read gauge transformation
    qlat::GaugeTransform qlat_gtinv;
    {
      qlat::GaugeTransform qlat_gt;
      qlat::dist_read_field(qlat_gt, gt_path);
      qlat::to_from_big_endian_64(qlat::get_data(qlat_gt)); 
      qlat::gt_inverse(qlat_gtinv, qlat_gt);
    }
    LatticeColourMatrix gt(grid);
    grid_convert(gt, qlat_gtinv);
    std::cout << "after gt" << std::endl;

    for(int t=0; t<T; ++t) {
      std::string path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/luchang/wall_s/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/strange ; t=" + std::to_string(t); // s quark wall source
      // string path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/luchang/gauge_transform_and_wall_l/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t); // l quark wall source

      std::string out_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/wall_s/" + std::to_string(traj) + "/" + std::to_string(t);
      // std::string out_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/wall_l/" + std::to_string(traj) + "/" + std::to_string(t);

      read_qlat_propagator(prop, path);

      prop = gt * prop;  // gt is the inverse of saved gauge transformation
      writeScidac_prop_d2f(prop, out_path);
    }
  }

  qlat::end();

  return 0;
}
