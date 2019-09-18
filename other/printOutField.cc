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

  LatticePGG ret(grid);
  readScidac(ret, "/projects/CSC249ADSE03/yidizhao/pGG_config/32ID/disc_2/pGG_disc.1250");
  std::cout << ret << std::endl;


  end();

  return 0;
}
