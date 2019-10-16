#include <headers/headers.h>
#include "../env.h"
#include "../util.h"


using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


std::vector<int> gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  Env env(gcoor, "24ID");

  LatticeKGG kgg(env.grid);
  readScidac(kgg, argv[1]);

  std::cout << "[1,1,1,1]" << std::endl;
  print_grid_field_site(kgg, {1,1,1,1});

  std::cout << "[1,2,3,4]" << std::endl;
  print_grid_field_site(kgg, {1,2,3,4});
  return 0;
}
