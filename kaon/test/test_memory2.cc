#include <map>
#include <iostream>
#include <fstream>
#include <string>
void print_memory() {
  using namespace std;

  map<string, double> memInfo; // item -> memory in GB
  ifstream f("/proc/meminfo");
  while(f) {
    string item, tmp;
    long long int size;

    f >> item >> size;
    getline(f, tmp); // move to the next line

    if(!item.empty()) {
      item.pop_back(); // remove last character
      memInfo[item] = double(size) / 1024 / 1024;
    }
  }
  // for(auto [s, i]: memInfo) std::cout << s << " " << i << " KB" << std::endl;

  std::cout << "MemTotal: " << memInfo["MemTotal"] << " GB " << std::endl;
  std::cout << "Non cache/buffer MemUsed: " << memInfo["MemTotal"] - memInfo["MemFree"] - memInfo["Buffers"] - memInfo["Cached"] << " GB " << std::endl;
  std::cout << "Buffers/Cached: " << memInfo["Buffers"] + memInfo["Cached"] << " GB " << std::endl;
  std::cout << "MemAvailable: " << memInfo["MemAvailable"] << " GB " << std::endl; // Memavailable is roughly memfree + buffers + cached,
}



#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;



int main(int argc, char* argv[])
{
  // std::cout << "Before Grid init" << std::endl;
  Grid_init(&argc, &argv);
  // std::cout << "after Grid init" << std::endl;
  // Coordinate mpi_coor = GridDefaultMpi();
  // std::cout << "after xxxxx" << std::endl;


  Coordinate gcoor({24, 24, 24, 64});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplexD::Nsimd()), GridDefaultMpi()); 

  grid->show_decomposition();
  
  print_memory();
  std::cout << "starting to create vector" << std::endl;
  std::cout << string(40, '=') << std::endl;
  vector<LatticePropagator> vec(64, grid);
  std::cout << string(40, '=') << std::endl;
  std::cout << "finished creating vector" << std::endl;
  print_memory();
  // LatticePropagator prop(grid);

  std::cout << "starting to create vector2" << std::endl;
  std::cout << string(40, '=') << std::endl;
  vector<LatticePropagator> vec2(64, grid);
  std::cout << string(40, '=') << std::endl;
  std::cout << "finished creating vector2" << std::endl;
  print_memory();


  return 0;
}
