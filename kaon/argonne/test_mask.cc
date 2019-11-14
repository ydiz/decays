#include <headers/headers.h>
#include "env.h"
#include "util.h"


using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


std::vector<int> gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();

  Env env(gcoor, "24ID");

  LatticeComplex mask(env.grid);

  int t_wall = 10;
  zero_mask(mask, t_wall);

  std::cout << mask << std::endl;

  return 0;
}
