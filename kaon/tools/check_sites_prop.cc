// #include <headers/headers.h>
// #include "../env.h"
// #include "../util.h"
#include "../kaon.h"


using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


std::vector<int> gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  Env env(gcoor, "24ID");

  LatticePropagator lat(env.grid);
  readScidac_prop_f2d(lat, argv[1]);

  std::cout << "[1,1,1,1]" << std::endl;
  print_grid_field_site(lat, {1,1,1,1});

  std::cout << "[1,2,3,4]" << std::endl;
  print_grid_field_site(lat, {1,2,3,4});
  return 0;
}
