#include "connected.h"
#include "disc.h"

// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});

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

  string ensemble = "24ID";
  // int t_min = 10; // should be consistent with calculating loop2
  int t_min = 20; // should be consistent with calculating loop2
  double Mpi = 0.13975;

  int traj_start = 1000, traj_end = 2290, traj_sep = 10; // need to have both loop and wall src propagators 
  // int traj_start = 2240, traj_end = 2290, traj_sep = 10; // need to have both loop and wall src propagators 
  std::vector<int> traj_skip {1020, 1060, 1100, 1340};
  std::cout << "traj_start: " << traj_start << std::endl;
  std::cout << "traj_end: " << traj_end << std::endl;
  std::cout << "traj_sep: " << traj_sep << std::endl;
  std::cout << "traj_skip: " << traj_skip << std::endl;
  // int traj_num = (traj_end - traj_start) / traj_sep + 1;
  // std::vector<int> trajs(traj_num);
  // for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * traj_sep;
  // // std::vector<int> trajs {1370};
  //
  // cout << "trajs: " << endl;
  // cout << trajs << endl;

  // for(int traj: trajs) {
  for(int traj = traj_start; traj<=traj_end; traj+=traj_sep) {
    if(std::find(traj_skip.begin(), traj_skip.end(), traj) != traj_skip.end()) continue;
    std::cout << "traj: " << traj << std::endl;

    LatticePGG disc(grid);
    calculate_disc(disc, "24ID", traj, t_min, Mpi);
    writeScidac(disc, "/projects/CSC249ADSE03/yidizhao/pGG_config/" + ensemble + "/disc/t-min="+ std::to_string(t_min) + "/pGG_disc." + std::to_string(traj));

    std::cout.flush();

    // cout << disc << endl;
  }


  end();

  return 0;
}
