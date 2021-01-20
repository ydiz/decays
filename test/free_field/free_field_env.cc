// Must add -DCUTH_FREE_FIELD to makefile, and not add -DBNL
#include "../../kaon/kaon.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;
int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  Env env("FreeField_8nt8");

  // std::vector<LatticePropagator> wl = env.get_wall('l');
  // std::cout << wl[5] << std::endl;

  // std::vector<LatticePropagator> ws = env.get_wall('s');
  // std::cout << ws[5] << std::endl;

  // LatticePropagator pl = env.get_point({1,2,3,4}, 'l');
  // std::cout << pl << std::endl;

  // LatticePropagator ps = env.get_point({1,2,3,4}, 's');
  // std::cout << ps << std::endl;

  LatticePropagator seq = env.get_sequential({1,2,3,4});
  std::cout << seq << std::endl;

  std::cout << "Finished!" << std::endl;
  Grid_finalize();
}
