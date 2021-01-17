#include <Grid/Grid.h>


using namespace Grid;
using namespace std;


template<class T>
void print_grid_field_site(const T &field, const std::vector<int> coor) {
  using namespace Grid;
  std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
  typename T::vector_object::scalar_object site;
  Coordinate c(coor);
  peekSite(site, field, c);
  std::cout << site << std::endl;
}


int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  string Umu_dir = "/hpcgpfs01/work/lqcd/etap/chulwoo/evec/24ID1Gev/configurations";

  // Coordinate fdims({24, 24, 24, 64});
  Coordinate fdims({4, 4, 4, 4});
  // vector<int> fdims = {24, 24, 24, 64};

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  GridParallelRNG pRNG(UGrid); pRNG.SeedFixedIntegers({1,2,3,4});

  int traj_start = 2000;
  int traj_end = 2000;
  int traj_sep = 10;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {

    LatticeGaugeField Umu(UGrid);
    // FieldMetaData header;
    // NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));
    Umu = 1.0;


    // generate a random gauge transformation field g
    LatticeColourMatrix g(UGrid); 
    LatticeColourMatrix tmp(UGrid); // random su(3) matrix
    SU<3>::GaussianFundamentalLieAlgebraMatrix(pRNG, tmp);
    autoView(g_v, g, AcceleratorWrite);
    autoView(tmp_v, tmp, AcceleratorRead);
    accelerator_for(ss, g.Grid()->oSites(), 1, {
      g_v[ss] = ProjectOnGroup(Exponentiate(tmp_v[ss], 1.0));
    });
    print_grid_field_site(g, {0,1,2,3});


    SU<3>::GaugeTransform(Umu, g); // g(x) * U(x) * adj(g(x+mu))
    print_grid_field_site(Umu, {0,1,2,3});

    Real plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
    std::cout << " Initial plaquette "<< plaq << std::endl;

    // std::cout << Umu << std::endl;

    ////////////////////////////////////////
    // Option 2. Calculate gauge transformation
    ////////////////////////////////////////
    // double alpha = 0.1;
    double alpha = 0.01;  // zyd: set this to a small value
    int coulomb_dir = Nd - 1;
    LatticeColourMatrix tmp_gt(UGrid);
    LatticeGaugeField Umu_coulomb = Umu;
    FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu, tmp_gt, alpha, 10000, 1.0e-12, 1.0e-12, true, coulomb_dir);    // Coulomb gauge fixing
    // FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu_coulomb, tmp_gt, alpha, 10000, 1.0e-12, 1.0e-12, true);  // Landau gauge fixing   

    std::cout << "after gauge-fixing" << std::endl;
    print_grid_field_site(tmp_gt, {0,1,2,3});
    print_grid_field_site(Umu_coulomb, {0,1,2,3});


    // ??? tmp_gt is different from adj(g); but they both transform Umu to original cold configuration
    std::cout << "transform Umu with tmp_gt" << std::endl;
    SU<3>::GaugeTransform(Umu, tmp_gt); // Original Umu, transformed by tmp_gt, should equal to Coulomb gauge Umu. // However, tmp_gt is not equal to g or adj(g)
    print_grid_field_site(Umu, {0,1,2,3});

    // std::cout << "transform Umu with tmp_gt" << std::endl;
    // LatticeColourMatrix g_adj = adj(g);
    // SU<3>::GaugeTransform(Umu, g_adj); // Original Umu, transformed by tmp_gt, should equal to Coulomb gauge Umu. // However, tmp_gt is not equal to g or adj(g)
    // print_grid_field_site(Umu, {0,1,2,3});

    plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
    std::cout << " final plaquette "<< plaq << std::endl;
  }


  std::cout << "FINISHED!" << std::endl;
  Grid_finalize();


}
