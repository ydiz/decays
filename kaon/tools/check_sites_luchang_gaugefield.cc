#include "../kaon.h"


using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;


template<typename T>
void cmdOptionToVector(const std::string &str,std::vector<T> &vec)
{
  vec.resize(0);
  std::stringstream ss(str);
  T i;
  while (ss >> i){
    vec.push_back(i);
    if(std::ispunct(ss.peek()))
      ss.ignore();
  }
  return;
}


int main(int argc, char* argv[])
{
  // Grid_init(&argc, &argv);
  zyd_init_Grid_Qlattice(argc, argv);

  Env env("24ID");

  LatticeGaugeField lat(env.grid);

  if(argc < 2) {
    std::cout << "argc must be larger than 1" << std::endl;
  }

  read_luchang_dist_gaugefield(lat, argv[1]);



  // // REMOVE ME // For test
  // Real plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(lat);
  // std::cout << "plaquette "<< plaq << std::endl;
  //
  // double alpha = 0.01;  // zyd: set this to a small value
  // int coulomb_dir = Nd - 1;
  // LatticeColourMatrix tmp_gt(lat.Grid());
  // FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(lat, tmp_gt, alpha, 10000, 1.0e-3, 1.0e-3, true, coulomb_dir);  


  if(argc == 2) {
    std::cout << "[1,1,1,1]" << std::endl;
    print_grid_field_site(lat, {1,1,1,1});

    std::cout << "[1,2,3,4]" << std::endl;
    print_grid_field_site(lat, {1,2,3,4});
  }
  else if(argc > 2) {
    for(int i=2; i<argc; ++i) {

      std::vector<int> coor;
      cmdOptionToVector(argv[i], coor);
      std::cout << coor << std::endl;
      print_grid_field_site(lat, coor);
    }
  }
  return 0;
}
