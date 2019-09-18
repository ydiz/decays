#include "connected.h"
#include "disc.h"

// std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> gcoor({24, 24, 24, 64});

namespace Grid {
namespace QCD {
void read_three_point_connected(LatticePGG &three_point, const std::string &ensemble, int traj) {
  if(ensemble == "32ID") {
    // std::string file = three_point_exact_path(traj); 
    std::string file = three_point_32ID(traj); 
    read_cheng_PGG(three_point, file);
  }
  else if(ensemble == "32IDF") {
    std::string file = three_point_path_32IDF(traj);
    read_luchang_PGG(three_point, file); // FIXME: change this after cheng generated his three point functions
  }
  else if(ensemble == "48I") {
    std::string file = three_point_48I(traj);
    read_luchang_PGG(three_point, file); // FIXME: change this after cheng generated his three point functions
  }
  else if(ensemble == "24ID") {
    std::string file = three_point_24ID(traj);
    read_cheng_PGG(three_point, file);
  }
  else assert(0);
}

}}


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
  int t_min_disc = 10; // should be consistent with calculating loop2
  double Mpi = 0.14;

  // int traj_start = 680, traj_num = 70;
  // std::vector<int> trajs(traj_num);
  // for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;
  std::vector<int> trajs {1370};

  cout << "trajs: " << endl;
  cout << trajs << endl;

  for(int traj: trajs) {

    LatticePGG disc(grid);
    calculate_disc(disc, ensemble, traj, t_min_disc, Mpi);

    LatticePGG conn(grid);
    read_three_point_connected(conn, ensemble, traj);
    conn = imag(conn); //  Cheng multiplied i to pion interpolator. 
    LatticeComplex tmp(grid);
    get_translational_factor(tmp, Mpi);
    conn = tmp * conn;

    LatticePGG threePoint(grid);
    threePoint = conn - disc; // connected and disconnected graph has a relative minus sign

    writeScidac(threePoint, "/home/yidizhao/cooley/pGG_config/" + ensemble + "/pGG." + std::to_string(traj));
    cout << threePoint << endl;
  }


  end();

  return 0;
}
