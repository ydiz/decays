// 3-dim grid is Okay to use; see  Grid/tests/core/Test_where.cc 

#include <Grid/Grid.h>
#include "kaon/kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  grid->show_decomposition();

  LatticeKGG lep(grid);
  double M_K = 0.5;
  std::vector<int> v = {1,2,6,5};
  int max_uv_sep = 3;
  EM_factor_half_lattice(lep, v, M_K, max_uv_sep);
  std::cout << lep << std::endl;


  Grid_finalize();
}


