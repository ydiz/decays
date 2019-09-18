#include "disc2.h"

std::vector<int> gcoor({32, 32, 32, 64});
// std::vector<int> gcoor({24, 24, 24, 64});

using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

	// int traj_start = 1180;
  double Mpi = 0.139474;
  int traj_start = std::stoi(GridCmdOptionPayload(argv,argv+argc,"--traj_start"));
  int traj_num = std::stoi(GridCmdOptionPayload(argv,argv+argc,"--traj_num"));
  std::vector<int> trajs(traj_num);
  for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;

  cout << "trajs: " << endl;
  cout << trajs << endl;

  for(int traj: trajs) {

    generatePointLoop("32ID", traj);

    // LatticePGG ret(grid);
    // three_point_disc2(ret, "32ID", traj, Mpi);
    // writeScidac(ret, "/projects/CSC249ADSE03/yidizhao/pGG_config/32ID/disc_2/pGG_disc." + std::to_string(traj));
  }


  end();

  return 0;
}
