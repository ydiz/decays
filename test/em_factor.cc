// 3-dim grid is Okay to use; see  Grid/tests/core/Test_where.cc 

#include <Grid/Grid.h>
#include "kaon/kaon.h"
#include "amplitude/form_factor.h"

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
  form_factor_integrand(lep, M_K);
  std::cout << lep << std::endl;


  Grid_finalize();
}


