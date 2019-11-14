#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main(int argc, char* argv[])
{
  std::cout << "Before Grid init" << std::endl;
  Grid_init(&argc, &argv);
  std::cout << "after Grid init" << std::endl;
  Coordinate mpi_coor = GridDefaultMpi();
  std::cout << "after xxxxx" << std::endl;


Coordinate gcoor({24, 24, 24, 64});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplexD::Nsimd()), mpi_coor); 
  LatticePropagator prop(grid);
  return 0;
}
