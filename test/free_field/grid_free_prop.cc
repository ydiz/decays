#include <Grid/Grid.h>
#include "kaon/kaon.h"

using namespace std;
using namespace Grid;


void print_free_field_prop(const LatticePropagator &prop) {
  LatticeSpinMatrix tmp = peekColour(prop, 0, 0);
  std::cout << tmp << std::endl;
}

int main(int argc, char **argv)
{
  Grid_init(&argc,&argv);


  Coordinate fdims({8,8,8,8});
  int Ls=16;
  double mass=0.04, M5=1.2;
  double time_boundary_phase = 113. * M_PI / 180.;  // 113 degrees
  // double mass=0.04, M5=1.0;  // M5 cannot be 1.0; result will be nan

  // GridCartesian         GRID(latt_size, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  // GridRedBlackCartesian RBGRID(&GRID);
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  LatticeGaugeFieldD Umu(UGrid); // Not used; but needed to create Ddwf
  Umu = 1.0;
  DomainWallFermionD Ddwf(Umu,*FGrid,*FrbGrid, *UGrid, *UrbGrid, mass, M5);


  // Point source propagator
  // Coordinate point_src({0,0,0,0});
  Coordinate point_src({1,2,3,4});  // When source is not 0, Grid free propagator code is wrong.

  LatticePropagator point_prop(UGrid);
  for(int s = 0; s < 4; ++s) {
    for(int c = 0; c < 3; ++c) {
      std::cout << GridLogMessage << "s: " << s << " c: " << c << std::endl;
      LatticeFermion src(UGrid), sol(UGrid);
      // Construct 4d source
      src = Zero();
      LatticeFermionSite site = Zero();
      site()(s)(c) = 1.;
      pokeSite(site, src, point_src);

      // double time_boundary_phase = 0.;  // 113 degrees
      std::vector<Complex> boundary(4, 1.0);
      boundary[3] = std::exp(Complex(0., time_boundary_phase)); // exp(i alpha)
      std::vector<double> twist(4, 0.);  // Always set twist to 0; only change boundary; 
      bool fiveD = false; //calculate 4d free propagator
      Ddwf.FreePropagator(src, sol, mass, boundary, twist, fiveD);

      FermToProp<WilsonImplR>(point_prop, sol, s, c);
    }
  }
  print_free_field_prop(point_prop);

  Grid_finalize();
}


