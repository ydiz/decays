#include <Grid/Grid.h>
#include "kaon/kaon.h"

using namespace std;
using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  grid->show_decomposition();

  double M_K = 0.5;
  int max_uv_sep = 3;

  // LatticeKGG lep(grid);
  // std::vector<int> v = {1,2,6,5};
  // EM_factor_half_lattice(lep, v, M_K, max_uv_sep);
  // std::cout << lep << std::endl;

  vector<vector<LatticeComplex>> Euv_fft_conj;
  Euv_fft_conj = calc_Euv_fft_conj(grid, M_K, max_uv_sep);
  std::cout << Euv_fft_conj[2][0] << std::endl;

  Grid_finalize();
}


