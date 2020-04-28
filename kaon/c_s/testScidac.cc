<<<<<<< HEAD
// #include <Grid/Grid.h>
=======
>>>>>>> 3c63df76e252c0fe6d63c6f2634fc315f64df691
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
<<<<<<< HEAD
  readScidac(lat, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/gauge_transform/2300");
=======
  readScidac(lat, "./2300");
>>>>>>> 3c63df76e252c0fe6d63c6f2634fc315f64df691
    std::cout << "finished" << std::endl;
  Grid_finalize();


  return 0;
}
