// #include <headers/headers.h>
// #include "../env.h"
// #include "../util.h"
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


std::vector<int> gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);

  Env env(gcoor, "24ID");

  LatticePropagator lat(env.grid);
  readScidac_prop_f2d(lat, argv[1]);


  if(argc < 2) {
    std::cout << "argc must be larger than 1" << std::endl;
  }
  else if(argc == 2) {
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
