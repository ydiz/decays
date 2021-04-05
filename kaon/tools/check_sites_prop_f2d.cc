// #include "../util.h"
#include "../kaon.h"


using namespace std;
using namespace Grid;


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
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({24,24,24,64}), GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  if(argc < 2) {
    std::cout << "argc must be larger than 1" << std::endl;
  }

  LatticePropagator lat(grid);
  readScidac_prop_f2d(lat, argv[1]);

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
