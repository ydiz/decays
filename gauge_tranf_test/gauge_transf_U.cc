// #include <Grid/Grid.h>
#include "../kaon/kaon.h"


using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;
using namespace std;


// template<class T>
// void print_grid_field_site(const T &field, const std::vector<int> &coor) {
//   using namespace Grid;
//   std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
//   typename T::vector_object::scalar_object site;
//   // Coordinate c(coor);
//   std::vector<int> c = coor;
//   peekSite(site, field, c);
//   std::cout << site << std::endl;
// }


int main(int argc, char **argv)
{
  // Grid_init(&argc, &argv);
  zyd_init_Grid_Qlattice(argc, argv);

  string Umu_dir = "/hpcgpfs01/work/lqcd/etap/chulwoo/evec/24ID1Gev/configurations";

  // Coordinate fdims({24, 24, 24, 64});
  // Coordinate fdims({4, 4, 4, 4});
  vector<int> fdims = {24, 24, 24, 64};

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(fdims, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());
  GridParallelRNG pRNG(grid); pRNG.SeedFixedIntegers({1,2,3,4});

  int traj_start = 2000;
  int traj_end = 2000;
  int traj_sep = 10;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {

    LatticeGaugeField Umu(grid);
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, Umu_dir + "/ckpoint_lat." + to_string(traj));
    // Umu = 1.0;


    // perform a random gauge transformation
    LatticeColourMatrix g(grid); 
    LatticeGaugeField Umu_new = Umu;
    SU<3>::RandomGaugeTransform(pRNG, Umu_new, g);

    string U_new_dir = "/sdcc/u/ydzhao/decays/gauge_tranf_test/saved_field/config/";
    NerscIO::writeConfiguration(Umu_new, U_new_dir + "/ckpoint_lat." + to_string(traj), 0, 0);

    // read luchang's g_to_Coul
    LatticeColourMatrix g_to_Coul(grid); // Luchang's g that sets Umu to Coulomb gauge
    string g_dir = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_luchang_props/wall_l/results/results=" + to_string(traj) + "/huge-data/gauge-transform";
    read_luchang_dist_gt(g_to_Coul, g_dir);

    // // Check Luchange's Coulomb gauge fixed configuration is g_to_Coul(x) * Umu(x) * g_to_Coul(x+mu)^dagger
    //
    // LatticeGaugeField Umu_Coul(grid);
    // string Umu_Coul_dir = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_luchang_props/wall_l/results/results=" + to_string(traj) + "/huge-data/gauge-fixed-gauge-field";
    // read_luchang_dist_gaugefield(Umu_Coul, Umu_Coul_dir);
    //
    // SU<3>::GaugeTransform(Umu, g_to_Coul);
    //
    // print_grid_field_site(Umu_Coul, {0,1,2,3});
    // print_grid_field_site(Umu, {0,1,2,3});
    // LatticeGaugeField tmp = Umu_Coul - Umu;
    // std::cout << norm2(tmp) << std::endl;


    LatticeColourMatrix g_to_Coul_new = g_to_Coul * adj(g); // Luchang's g that sets Umu to Coulomb gauge
    string g_new_dir = "/sdcc/u/ydzhao/decays/gauge_tranf_test/saved_field/gauge_transformation/";
    writeScidac(g_to_Coul_new, g_new_dir + "/gt." + to_string(traj));

    // // Check that new gauge transformation g_to_Coul_new indeed transforms Umu_new to Coulomb gauge
    LatticeGaugeField Umu_Coul(grid);
    string Umu_Coul_dir = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID_luchang_props/wall_l/results/results=" + to_string(traj) + "/huge-data/gauge-fixed-gauge-field";
    read_luchang_dist_gaugefield(Umu_Coul, Umu_Coul_dir);

    SU<3>::GaugeTransform(Umu_new, g_to_Coul_new);

    LatticeGaugeField tmp = Umu_Coul - Umu_new;
    std::cout << norm2(tmp) << std::endl;

  }


  std::cout << "FINISHED!" << std::endl;
  Grid_finalize();


}
