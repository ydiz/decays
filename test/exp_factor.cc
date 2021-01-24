#include <Grid/Grid.h>
#include "kaon/kaon.h"

using namespace std;
using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  grid->show_decomposition();

  int tK = 7;
  double M_K = 0.5;

  LatticeComplex exp_factor = exp_v0_tK(grid, tK, M_K);
  std::cout << exp_factor << std::endl;

  Grid_finalize();
}


