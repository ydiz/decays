#include "../kaon.h"

// using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

std::vector<int> gcoor({24, 24, 24, 64});



int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  LatticeColourMatrix lat(grid);
  lat = 1.;
  //   std::cout << "writing" << std::endl;
  // writeScidac(lat, "./xx");
    std::cout << "reading" << std::endl;
  // readScidac(lat, "./xx");
  readScidac(lat, "./2300");
    std::cout << "finished" << std::endl;
  Grid_finalize();


  return 0;
}
